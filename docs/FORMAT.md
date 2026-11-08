# Data and output formats

## Object payloads

The database and optional external-query file share this text format:

```text
N M
metric_type_0 ... metric_type_M-1
scale_0 ... scale_M-1
object_0_component_0 ... object_0_component_M-1
...
```

Each object occupies one line. Components are whitespace-separated. Numeric
vectors use comma-separated values without spaces. All records must have the
same component dimensions. Each scale is a finite positive divisor applied to
the corresponding raw distance. External-query metadata must match the database.
Use finite numeric inputs and nonempty strings; whitespace cannot occur inside
string or set fields.

| Type | Representation | Distance |
|---|---|---|
| 1 | Scalar or comma-separated vector | L1 |
| 2 | Comma-separated vector | Euclidean L2, not squared L2 |
| 5 | Nonzero comma-separated vector | Angular distance in degrees |
| 6 | Nonempty UTF-8 byte string without whitespace | Byte-wise Levenshtein edit distance |
| 7 | `token_a|token_b|token_c` | Jaccard distance on token sets |

The example generator uses types 2, 1, 6 and 7. Type 5 uses the historical
`acos(cosine) * 180 / 3.1415926` convention; choose the scale in the same units.
Legacy aliases 0 and 4 follow L1; type 3 is the historical date/time encoding.
They are not needed for the example. To add a distance, extend `dist` in
`src/metric.cpp` and provide its corresponding input representation.

## Query protocol

```text
weight_0 ... weight_M-1
number_of_k_values
k_0 ... k_last
0
number_of_queries
query_id_0 ... query_id_last
```

The `0` line is the number of legacy range-query thresholds; keep it zero for
k-NN search. IDs are zero-based. With external queries they index the query
payload file; otherwise they index the database. Example:

```text
0.25 0.25 0.25 0.25
1
10
0
3
0 1 2
```

Use finite nonnegative weights, with at least one positive weight. The research
implementation uses small activity thresholds (approximately 1e-6 in EF search),
so use exact zero for disabled components and ordinary-scale positive weights.
Weights are not automatically normalized; use sum-one weights for consistency.
The requested k should not exceed the database size. Reuse the index when changing
weights or k.

## Answers

Each query emits:

```text
k <largest_returned_distance> <count> <database_id> ...
D <distance> ...
```

`D` appears only with `OUT_DIST=1`. IDs always refer to database rows, including
external-query mode. Each requested k is followed by:

```text
t <mean_time_ms>
dc <mean_component_distance_evaluations>
```

The exact-scan routine may emit neighbors in descending distance order; compare
sets or sort by distance rather than assuming the same order as graph search.
Ties at rank k can produce different valid result sets.

## Indexes

`build` emits a legacy text `.ind` and a native binary `.sjg` index. Use `.sjg`
for the documented workflow. It stores adjacency, entry information and per-edge
component distances. Object payloads remain in the database file. The binary
format uses the implementation's native fixed-width layout; use the same
platform/endianness and rebuild indexes when changing the payload or metrics.
