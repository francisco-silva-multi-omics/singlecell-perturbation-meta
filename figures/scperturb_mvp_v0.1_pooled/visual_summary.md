# scPerturb MVP Visual Summary

## Dataset Summary

| Dataset | Cells | Genes | Contrasts | Pass | Warn | Controls |
|---|---:|---:|---:|---:|---:|---:|
| adamson | 65,257 | 32,738 | 96 | 0 | 96 | n/a |
| dixit | 28,034 | 23,111 | 26 | 0 | 26 | n/a |
| norman | 111,445 | 33,694 | 236 | 0 | 236 | n/a |
| replogle | 310,385 | 8,563 | 1,832 | 1,832 | 0 | n/a |
| sciplex3 | 762,795 | 110,983 | 2,423 | 2,423 | 0 | n/a |

## Effect Magnitude

| Dataset | Contrasts | Median Max Abs Effect | P90 Max Abs Effect | Median P95 Abs Effect |
|---|---:|---:|---:|---:|
| adamson | 96 | 1.887 | 2.715 | 0.216 |
| dixit | 26 | 1.103 | 1.924 | 0.175 |
| norman | 236 | 3.056 | 4.360 | 0.333 |
| replogle | 1832 | 1.924 | 2.583 | 0.374 |
| sciplex3 | 2423 | 3.445 | 4.809 | 0.838 |

## Recurrent Top-Effect Genes

| Gene | Datasets | Top-Gene Hits | Max Abs Effect |
|---|---:|---:|---:|
| TMEM158 | 5 | 805 | 4.030 |
| AREG | 5 | 489 | 3.800 |
| IER3 | 5 | 473 | 4.431 |
| EGR1 | 5 | 445 | 4.174 |
| TUBA1A | 5 | 391 | 3.366 |
| TAC3 | 5 | 377 | 3.897 |
| TSPAN32 | 5 | 324 | 2.930 |
| NTRK1 | 5 | 319 | 3.376 |
| CAMK2N1 | 5 | 278 | 4.043 |
| HOMER3 | 5 | 266 | 2.838 |
| TRANK1 | 5 | 259 | 3.869 |
| ZNF467 | 5 | 252 | 2.910 |
| HSPA1A | 5 | 211 | 3.387 |
| S100A4 | 5 | 200 | 3.507 |
| ARHGEF6 | 5 | 199 | 4.170 |

## Figures

- overview: `figures\scperturb_mvp_v0.1_pooled\mvp_overview.svg`
- effect_magnitude: `figures\scperturb_mvp_v0.1_pooled\effect_magnitude.svg`
- top_gene_recurrence: `figures\scperturb_mvp_v0.1_pooled\top_gene_recurrence.svg`
