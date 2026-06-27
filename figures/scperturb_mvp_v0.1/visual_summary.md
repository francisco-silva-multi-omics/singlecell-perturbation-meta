# scPerturb MVP Visual Summary

## Dataset Summary

| Dataset | Cells | Genes | Contrasts | Controls |
|---|---:|---:|---:|---:|
| adamson | 65,257 | 32,738 | 96 | 1 |
| dixit | 28,034 | 23,111 | 26 | 4 |
| norman | 111,445 | 33,694 | 286 | 4 |
| replogle | 310,385 | 8,563 | 1,922 | 97 |
| sciplex3 | 762,795 | 110,983 | 188 | 1 |

## Effect Magnitude

| Dataset | Contrasts | Median Max Abs Effect | P90 Max Abs Effect | Median P95 Abs Effect |
|---|---:|---:|---:|---:|
| adamson | 96 | 1.887 | 2.715 | 0.216 |
| dixit | 26 | 1.103 | 1.924 | 0.175 |
| norman | 286 | 3.040 | 4.342 | 0.342 |
| replogle | 1922 | 1.937 | 2.575 | 0.378 |
| sciplex3 | 188 | 1.605 | 3.356 | 0.199 |

## Recurrent Top-Effect Genes

| Gene | Datasets | Top-Gene Hits | Max Abs Effect |
|---|---:|---:|---:|
| TMEM158 | 5 | 826 | 3.240 |
| AREG | 5 | 496 | 3.468 |
| EGR1 | 5 | 449 | 2.383 |
| TUBA1A | 5 | 399 | 3.365 |
| TAC3 | 5 | 343 | 2.094 |
| TSPAN32 | 5 | 334 | 2.044 |
| HSPA1A | 5 | 206 | 2.591 |
| ARHGEF6 | 5 | 179 | 2.665 |
| CAB39L | 5 | 156 | 1.604 |
| HIST1H3H | 5 | 140 | 1.728 |
| PRSS57 | 5 | 138 | 1.703 |
| RNF208 | 5 | 95 | 2.117 |
| VGF | 4 | 769 | 2.509 |
| CRYM | 4 | 625 | 2.859 |
| ALAS2 | 4 | 600 | 2.779 |

## Figures

- overview: `figures\scperturb_mvp_v0.1\mvp_overview.svg`
- effect_magnitude: `figures\scperturb_mvp_v0.1\effect_magnitude.svg`
- top_gene_recurrence: `figures\scperturb_mvp_v0.1\top_gene_recurrence.svg`
