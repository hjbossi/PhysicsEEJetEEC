# Data for the publication figure

| file | description |
|---|---|
| `points_Sept15th.txt` | Unfolded ALEPH EEC vs `z = (1 - cos θ)/2`, with statistical and systematic uncertainties (delivered 2026-09-15; supersedes `points_Oct30th.txt`, which had the same central values and systematics but older statistical uncertainties). |
| `5theory_bands_with_evo_omega1g_zstar004.txt` | NNLL_col + NNLO_FO + NNNNLL_b2b prediction: central, lower and upper tables as Mathematica lists (`expfulltable[...] = {{z, EEC}, ...}`). |

## `points_Sept15th.txt` format

One line per bin:

```
index: N { z_low , EEC , sys: S stat: T }
```

* `N` – bin number on the analysis grid (see below); the file contains bins 17–183.
* `z_low` – **lower edge** of bin `N`. This is **not** the bin centre. The upper
  edge of bin `N` is the lower edge of bin `N+1`; for the last line it must be
  taken from the grid.
* `EEC` – bin content (normalised, divided by the bin width in `z`).
* `S`, `T` – systematic and statistical uncertainty on `EEC`.

Every line, including the last one, is a real bin with a real measurement.

## Analysis binning

200 bins in `z`: 100 logarithmically spaced bins from `1e-6` to `0.5`, mirrored
about `z = 0.5` so that the upper half is logarithmic in `1 - z`. In Python:

```python
edges = calcBinEdge(1e-6, 0.5, 100)   # 201 edges, see script/EEC_Plot_Pub_ALEPH.ipynb
```

The `z_low` column agrees with `edges[N]` to better than `1e-6` (relative).
The plotting notebook rebuilds the grid, checks the file against it, and uses
arithmetic bin midpoints as plotting positions. The theory tables are
tabulated at the same `z` values as the lower bin edges; the prediction is a
continuous curve and is drawn at those `z`.
