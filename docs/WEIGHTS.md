# Experimental query weights

[Back to README](../README.md) · [Download full-precision CSV](weights.csv)

The six datasets use seven fixed query-weight profiles: one uniform profile,
three dense profiles, and three skewed profiles. The same SJG graph is reused
across profiles; SJG construction does not use these query weights.

## Profiles

| Profile | Definition |
|---|---|
| `uniform` | Equal weights; the recorded protocol uses all ones |
| `dense_r0`–`dense_r2` | Three fixed draws from a symmetric Dirichlet distribution with concentration 1 |
| `skewed_r0`–`skewed_r2` | Three fixed draws with concentration 0.2 |

The nonuniform vectors sum to one. Uniform weights are reproduced as recorded,
rather than rescaled to sum to one. The weighted distance is
`sum_i weight_i * raw_distance_i / scale_i`; scales belong to the dataset metadata.
The `r0`–`r2` suffixes identify weight profiles, not timing repetitions.

Tables display six significant digits for readability. **Use the CSV values for
reproduction**, especially for very small nonzero weights. The CSV uses one row
per component, with a one-based component index. Vector positions follow the
component order listed below. The numeric values have been checked against the
saved experimental query protocols.

## Dataset profiles

### Corel

Dataset key: `corel`. Component order:

1. Color histogram (32D)
2. Color layout histogram (32D)
3. Color moments (9D)
4. Co-occurrence texture (16D)

| Profile | w1 | w2 | w3 | w4 |
|---|---:|---:|---:|---:|
| `uniform` | 1 | 1 | 1 | 1 |
| `dense_r0` | 0.667618 | 0.199062 | 0.0998866 | 0.0334334 |
| `dense_r1` | 0.109461 | 0.0615008 | 0.824436 | 0.00460208 |
| `dense_r2` | 0.0402514 | 0.260209 | 0.0199657 | 0.679574 |
| `skewed_r0` | 0.193985 | 0.000169038 | 0.520867 | 0.28498 |
| `skewed_r1` | 0.0517263 | 0.00234894 | 0.884373 | 0.0615519 |
| `skewed_r2` | 0.157524 | 0.0266601 | 0.000196881 | 0.815619 |

### NUS-WIDE

Dataset key: `nuswide`. Component order:

1. Color histogram (64D)
2. Color correlogram (144D)
3. Edge direction histogram (73D)
4. Wavelet texture (128D)
5. Block-wise color moments (225D)

| Profile | w1 | w2 | w3 | w4 | w5 |
|---|---:|---:|---:|---:|---:|
| `uniform` | 1 | 1 | 1 | 1 | 1 |
| `dense_r0` | 0.171913 | 0.289808 | 0.427469 | 0.0728741 | 0.0379357 |
| `dense_r1` | 0.334689 | 0.317815 | 0.17518 | 0.0960724 | 0.0762443 |
| `dense_r2` | 0.189413 | 0.0366524 | 0.0418653 | 0.414247 | 0.317822 |
| `skewed_r0` | 0.00583508 | 0.0010953 | 0.0162176 | 0.0138326 | 0.963019 |
| `skewed_r1` | 0.0116044 | 0.795451 | 0.156488 | 0.00346439 | 0.0329919 |
| `skewed_r2` | 0.0879123 | 0.301001 | 0.162894 | 0.140646 | 0.307547 |

### Tiny ImageNet

Dataset key: `imagenet`. Component order:

1. HSV histogram (256D)
2. HSV color moments (9D)
3. Gabor (32D)
4. HOG (324D)
5. ResNet-50 (2048D)

| Profile | w1 | w2 | w3 | w4 | w5 |
|---|---:|---:|---:|---:|---:|
| `uniform` | 1 | 1 | 1 | 1 | 1 |
| `dense_r0` | 0.570348 | 0.08977 | 0.0734327 | 0.149658 | 0.116792 |
| `dense_r1` | 0.0928703 | 0.351963 | 0.195266 | 0.193466 | 0.166435 |
| `dense_r2` | 0.139037 | 0.156679 | 0.203963 | 0.243299 | 0.257023 |
| `skewed_r0` | 0.00135455 | 1.0683e-06 | 0.136926 | 0.845709 | 0.0160091 |
| `skewed_r1` | 0.0660313 | 1.27858e-05 | 0.933236 | 0.000581186 | 0.00013865 |
| `skewed_r2` | 0.0191367 | 0.942557 | 0.0341203 | 9.22035e-05 | 0.00409347 |

### COCO

Dataset key: `coco`. Component order:

1. HSV histogram
2. Block-wise color moments
3. LBP
4. Edge direction histogram
5. CLIP image
6. ResNet-50 image
7. MPNet text
8. CLIP text

| Profile | w1 | w2 | w3 | w4 | w5 | w6 | w7 | w8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `uniform` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| `dense_r0` | 0.182187 | 0.182701 | 0.0638808 | 0.0568084 | 0.209603 | 0.0352806 | 0.0576378 | 0.211901 |
| `dense_r1` | 0.276807 | 0.0510764 | 0.0411062 | 0.115209 | 0.00620462 | 0.0772822 | 0.0856604 | 0.346654 |
| `dense_r2` | 0.137427 | 0.0967999 | 0.0147801 | 0.147626 | 0.174279 | 0.0604113 | 0.115891 | 0.252786 |
| `skewed_r0` | 0.184977 | 2.10837e-08 | 0.152797 | 0.0739351 | 0.168562 | 0.014515 | 0.0863552 | 0.31886 |
| `skewed_r1` | 0.00247298 | 0.132856 | 0.00591833 | 0.266603 | 0.213425 | 0.307914 | 0.0563595 | 0.0144513 |
| `skewed_r2` | 0.00604827 | 0.00942694 | 0.831005 | 0.0330324 | 2.03509e-06 | 0.117575 | 0.00167569 | 0.00123475 |

### Food

Dataset key: `foodhetero`. Component order:

1. Energy
2. Fat
3. Sugar
4. Salt
5. Protein
6. Food category
7. Product-name embedding

| Profile | w1 | w2 | w3 | w4 | w5 | w6 | w7 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `uniform` | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| `dense_r0` | 0.543575 | 0.0209462 | 0.0463235 | 0.130129 | 0.0642647 | 0.00919985 | 0.185562 |
| `dense_r1` | 0.0588659 | 0.10638 | 0.212741 | 0.14019 | 0.292299 | 0.0176913 | 0.171833 |
| `dense_r2` | 0.194106 | 0.105269 | 0.123996 | 0.168169 | 0.228136 | 0.0819135 | 0.0984106 |
| `skewed_r0` | 0.0220814 | 0.877042 | 0.0078451 | 0.0864767 | 2.40835e-07 | 5.18206e-05 | 0.00650307 |
| `skewed_r1` | 0.0948079 | 3.46125e-06 | 0.186568 | 0.0118427 | 0.7066 | 0.000177178 | 5.10691e-07 |
| `skewed_r2` | 0.000281999 | 1.53502e-05 | 0.538547 | 0.244052 | 4.77689e-05 | 0.217033 | 2.2079e-05 |

### Airbnb

Dataset key: `rental1m`. Component order:

1. Latitude/longitude
2. Listing name
3. Price
4. Guest capacity
5. Amenities

| Profile | w1 | w2 | w3 | w4 | w5 |
|---|---:|---:|---:|---:|---:|
| `uniform` | 1 | 1 | 1 | 1 | 1 |
| `dense_r0` | 0.0503368 | 0.313732 | 0.367157 | 0.137856 | 0.130919 |
| `dense_r1` | 0.198304 | 0.0245027 | 0.579334 | 0.197793 | 6.61536e-05 |
| `dense_r2` | 0.353889 | 0.149405 | 0.299641 | 0.0811322 | 0.115932 |
| `skewed_r0` | 0.166126 | 0.0114199 | 0.821764 | 7.79467e-07 | 0.00069 |
| `skewed_r1` | 9.77459e-05 | 0.0002744 | 0.00342546 | 0.830312 | 0.165891 |
| `skewed_r2` | 0.613858 | 0.345952 | 4.46112e-11 | 0.00405356 | 0.0361364 |

## Copy a vector into a query file

This prints a full-precision weight line for NUS-WIDE, using only the Python
standard library. Replace the first line of your query protocol with the output.

```bash
python3 - <<'PY'
import csv
with open('docs/weights.csv', newline='') as f:
    rows = [r for r in csv.DictReader(f)
            if r['dataset'] == 'nuswide' and r['profile'] == 'dense_r0']
rows.sort(key=lambda r: int(r['component_index']))
print(' '.join(r['weight'] for r in rows))
PY
```

See [the query format](FORMAT.md#query-protocol) for k and query IDs.
The synthetic example has its own four-component schema; benchmark weights
must be paired with the matching benchmark component order.
