from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
import os

import h5py
import numpy as np


class CurationError(RuntimeError):
    """Raised when an input H5AD violates an expected curation contract."""


def decode_array(group: h5py.Group, name: str) -> np.ndarray:
    obj = group[name]
    if isinstance(obj, h5py.Dataset):
        values = obj[()]
        if values.dtype.kind in "SO":
            return np.asarray(
                [value.decode() if isinstance(value, bytes) else str(value) for value in values],
                dtype=object,
            )
        return values
    if obj.attrs.get("encoding-type") != "categorical":
        raise CurationError(f"Unsupported encoded column: {name}")
    categories = decode_array(obj, "categories")
    codes = obj["codes"][()]
    return np.asarray(
        [categories[code] if code >= 0 else None for code in codes], dtype=object
    )


def dataframe_index(group: h5py.Group) -> np.ndarray:
    index_name = group.attrs["_index"]
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    return decode_array(group, index_name)


def _dataset_options(source: h5py.Dataset, *, chunks=None) -> dict:
    options = {}
    if source.compression:
        options["compression"] = source.compression
        options["compression_opts"] = source.compression_opts
    if source.shuffle:
        options["shuffle"] = True
    if source.fletcher32:
        options["fletcher32"] = True
    if chunks is not None:
        options["chunks"] = chunks
    return options


def _copy_attrs(source, destination) -> None:
    for key, value in source.attrs.items():
        destination.attrs[key] = value


def _write_categorical(parent: h5py.Group, name: str, values: np.ndarray) -> None:
    if name in parent:
        del parent[name]
    group = parent.create_group(name)
    group.attrs["encoding-type"] = "categorical"
    group.attrs["encoding-version"] = "0.2.0"
    group.attrs["ordered"] = False

    nonmissing = sorted({str(value) for value in values if value is not None})
    lookup = {value: index for index, value in enumerate(nonmissing)}
    codes = np.asarray(
        [-1 if value is None else lookup[str(value)] for value in values], dtype=np.int32
    )
    categories = group.create_dataset(
        "categories", data=np.asarray(nonmissing, dtype=h5py.string_dtype("utf-8"))
    )
    categories.attrs["encoding-type"] = "string-array"
    categories.attrs["encoding-version"] = "0.2.0"
    code_dataset = group.create_dataset("codes", data=codes, compression="gzip")
    code_dataset.attrs["encoding-type"] = "array"
    code_dataset.attrs["encoding-version"] = "0.2.0"


def _write_array(parent: h5py.Group, name: str, values: np.ndarray) -> None:
    if name in parent:
        del parent[name]
    array = np.asarray(values)
    if array.dtype.kind in "OSU":
        dataset = parent.create_dataset(
            name,
            data=np.asarray([str(value) for value in array], dtype=h5py.string_dtype("utf-8")),
        )
        dataset.attrs["encoding-type"] = "string-array"
    else:
        dataset = parent.create_dataset(name, data=array, compression="gzip")
        dataset.attrs["encoding-type"] = "array"
    dataset.attrs["encoding-version"] = "0.2.0"


def _copy_filtered_dataframe(
    source: h5py.Group,
    destination_parent: h5py.Group,
    name: str,
    keep: np.ndarray,
    updates: dict[str, np.ndarray],
) -> None:
    destination = destination_parent.create_group(name)
    _copy_attrs(source, destination)
    row_count = len(keep)
    for column_name, obj in source.items():
        if isinstance(obj, h5py.Dataset):
            if obj.ndim == 1 and len(obj) == row_count:
                values = obj[()][keep]
                dataset = destination.create_dataset(
                    column_name,
                    data=values,
                    **_dataset_options(obj, chunks=True if obj.chunks else None),
                )
                _copy_attrs(obj, dataset)
            else:
                source.copy(column_name, destination)
        elif obj.attrs.get("encoding-type") == "categorical":
            group = destination.create_group(column_name)
            _copy_attrs(obj, group)
            obj.copy("categories", group)
            codes = obj["codes"][()][keep]
            code_source = obj["codes"]
            code_dataset = group.create_dataset(
                "codes",
                data=codes,
                **_dataset_options(code_source, chunks=True if code_source.chunks else None),
            )
            _copy_attrs(code_source, code_dataset)
        else:
            source.copy(column_name, destination)

    for column_name, values in updates.items():
        selected = np.asarray(values, dtype=object if np.asarray(values).dtype.kind in "OSU" else None)[keep]
        if selected.dtype.kind in "OSU":
            _write_categorical(destination, column_name, selected)
        else:
            _write_array(destination, column_name, selected)

    index_name = destination.attrs["_index"]
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    existing_order = [
        value.decode() if isinstance(value, bytes) else str(value)
        for value in destination.attrs.get("column-order", [])
    ]
    for column_name in updates:
        if column_name != index_name and column_name not in existing_order:
            existing_order.append(column_name)
    destination.attrs["column-order"] = np.asarray(
        existing_order, dtype=h5py.string_dtype("utf-8")
    )


def _contiguous_runs(indices: np.ndarray) -> Iterable[tuple[int, int]]:
    if not len(indices):
        return
    start = previous = int(indices[0])
    for value in indices[1:]:
        value = int(value)
        if value != previous + 1:
            yield start, previous + 1
            start = value
        previous = value
    yield start, previous + 1


def _create_sparse_vector(destination, name, source, size):
    chunks = source.chunks
    if chunks:
        chunks = (min(chunks[0], max(1, size)),)
    dataset = destination.create_dataset(
        name,
        shape=(size,),
        maxshape=(None,),
        dtype=source.dtype,
        **_dataset_options(source, chunks=chunks or True),
    )
    _copy_attrs(source, dataset)
    return dataset


def _filter_csr(source: h5py.Group, destination: h5py.Group, keep: np.ndarray) -> None:
    indptr = source["indptr"][()]
    selected_rows = np.flatnonzero(keep)
    lengths = np.diff(indptr)[selected_rows]
    new_indptr = np.empty(len(selected_rows) + 1, dtype=indptr.dtype)
    new_indptr[0] = 0
    np.cumsum(lengths, out=new_indptr[1:])
    new_nnz = int(new_indptr[-1])
    output_data = _create_sparse_vector(destination, "data", source["data"], new_nnz)
    output_indices = _create_sparse_vector(destination, "indices", source["indices"], new_nnz)

    write_at = 0
    for first_row, last_row in _contiguous_runs(selected_rows):
        first = int(indptr[first_row])
        last = int(indptr[last_row])
        count = last - first
        output_data[write_at : write_at + count] = source["data"][first:last]
        output_indices[write_at : write_at + count] = source["indices"][first:last]
        write_at += count

    indptr_dataset = destination.create_dataset("indptr", data=new_indptr)
    _copy_attrs(source["indptr"], indptr_dataset)


def _filter_csc(source: h5py.Group, destination: h5py.Group, keep: np.ndarray) -> None:
    indptr = source["indptr"][()]
    old_nnz = len(source["data"])
    output_data = _create_sparse_vector(destination, "data", source["data"], old_nnz)
    output_indices = _create_sparse_vector(destination, "indices", source["indices"], old_nnz)
    row_map = np.full(len(keep), -1, dtype=source["indices"].dtype)
    row_map[keep] = np.arange(int(keep.sum()), dtype=row_map.dtype)
    new_indptr = np.empty(len(indptr), dtype=indptr.dtype)
    new_indptr[0] = 0

    write_at = 0
    for column in range(len(indptr) - 1):
        first, last = int(indptr[column]), int(indptr[column + 1])
        rows = source["indices"][first:last]
        selected = keep[rows]
        count = int(selected.sum())
        if count:
            output_indices[write_at : write_at + count] = row_map[rows[selected]]
            output_data[write_at : write_at + count] = source["data"][first:last][selected]
            write_at += count
        new_indptr[column + 1] = write_at

    output_data.resize((write_at,))
    output_indices.resize((write_at,))
    indptr_dataset = destination.create_dataset("indptr", data=new_indptr)
    _copy_attrs(source["indptr"], indptr_dataset)


def _copy_filtered_matrix(
    source: h5py.Dataset | h5py.Group,
    destination_parent: h5py.Group,
    keep: np.ndarray,
) -> None:
    if isinstance(source, h5py.Dataset):
        chunks = source.chunks
        if chunks:
            chunks = (min(chunks[0], max(1, int(keep.sum()))), chunks[1])
        destination = destination_parent.create_dataset(
            "X",
            shape=(int(keep.sum()), source.shape[1]),
            dtype=source.dtype,
            **_dataset_options(source, chunks=chunks),
        )
        _copy_attrs(source, destination)
        write_at = 0
        for first, last in _contiguous_runs(np.flatnonzero(keep)):
            count = last - first
            destination[write_at : write_at + count] = source[first:last]
            write_at += count
        return

    destination = destination_parent.create_group("X")
    _copy_attrs(source, destination)
    encoding = source.attrs.get("encoding-type")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding == "csr_matrix":
        _filter_csr(source, destination, keep)
    elif encoding == "csc_matrix":
        _filter_csc(source, destination, keep)
    else:
        raise CurationError(f"Unsupported X encoding: {encoding}")
    shape = np.asarray(source.attrs["shape"]).copy()
    shape[0] = int(keep.sum())
    destination.attrs["shape"] = shape


def filter_h5ad_rows(
    input_path: Path,
    output_path: Path,
    keep: np.ndarray,
    *,
    obs_updates: dict[str, np.ndarray] | None = None,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".partial")
    if output_path.exists():
        raise CurationError(f"Output already exists: {output_path}")
    if temporary.exists():
        temporary.unlink()

    try:
        with h5py.File(input_path, "r") as source, h5py.File(temporary, "w") as destination:
            keep = np.asarray(keep, dtype=bool)
            if len(keep) != len(dataframe_index(source["obs"])):
                raise CurationError("Cell mask length does not match obs")
            _copy_attrs(source, destination)
            for key in source.keys():
                if key not in {"X", "obs"}:
                    source.copy(key, destination)
            _copy_filtered_matrix(source["X"], destination, keep)
            _copy_filtered_dataframe(
                source["obs"], destination, "obs", keep, obs_updates or {}
            )
        os.replace(temporary, output_path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
