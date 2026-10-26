# Representative MLIPPER Datasets

These fixtures cover different tip and site sizes without adding the complete
ROADIES experiment to Git history.

## Small-tip workflow (ROADIES)

| Case | Source | Reference tips | Query tips | Total tips | Sites |
| --- | --- | ---: | ---: | ---: | ---: |
| `small_short` | iteration 4, gene 1561 | 22 | 29 | 51 | 543 |
| `large_short` | iteration 1, gene 640 | 240 | 227 | 467 | 543 |
| `medium` | iteration 5, gene 6250 | 85 | 95 | 180 | 632 |
| `small_long` | iteration 4, gene 2558 | 23 | 26 | 49 | 772 |
| `large_long` | iteration 4, gene 1152 | 234 | 241 | 475 | 760 |

Every case contains only the inputs needed to run the workflow:
`reference.fa`, `query.fa`, `backbone.nwk`, and `model.bestModel`.

## Divide-and-conquer workflow (DIPPER)

| Case | Tips | Sites | Compressed | Extracted |
| --- | ---: | ---: | ---: | ---: |
| `10k` | 10,000 | 16,092 | 2.35 MB | 164 MB |
| `20k` | 20,000 | 20,812 | 6.49 MB | 423 MB |

These are the original aligned DIPPER AliSim fixtures, stored as gzip archives.
Each directory also contains the original simulation truth tree and a neutral
GTR+FC+G4 starting model for MLIPPER. Generated DIPPER and MLIPPER regression
outputs are not stored here. The model is not presented as the unknown original
AliSim simulation model; MLIPPER runs should pass `--empirical-freqs`.

Extract one or both archives into local output directories. These directories
remain Git-ignored; only `data/small_tip` and `data/divide_and_conquer` are
tracked:

```bash
mkdir -p data/dipper_10k data/dipper_20k
gzip -dc data/divide_and_conquer/10k/alignment.fa.gz \
  > data/dipper_10k/alignment.fa
gzip -dc data/divide_and_conquer/20k/alignment.fa.gz \
  > data/dipper_20k/alignment.fa
```

The DIPPER fixtures exercise MLIPPER's D&C input path; they are not ROADIES
placement fixtures and do not require a backbone tree.
