from __future__ import annotations

import json
import time
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


USER_AGENT = "singlecell-perturbation-meta/0.1 (+https://github.com/francisco-silva-multi-omics/singlecell-perturbation-meta)"


def open_with_retries(
    request: Request,
    *,
    timeout: float = 60,
    retries: int = 3,
):
    for attempt in range(retries + 1):
        try:
            return urlopen(request, timeout=timeout)
        except HTTPError as exc:
            retryable = exc.code == 429 or 500 <= exc.code < 600
            if not retryable or attempt == retries:
                raise
        except (TimeoutError, URLError):
            if attempt == retries:
                raise
        time.sleep(2**attempt)
    raise RuntimeError("unreachable")


def get_json(url: str, *, timeout: float = 60, retries: int = 3) -> dict[str, Any]:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    with open_with_retries(request, timeout=timeout, retries=retries) as response:
        return json.load(response)

