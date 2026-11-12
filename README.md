# SJG: Skyline Joint Graph

A compact C++17 implementation of SJG for approximate nearest-neighbor search
in multi-metric spaces. One graph serves queries with different nonnegative
component weights.

SJG builds component-wise nearest-neighbor candidates using NN-Descent, merges
forward and reverse neighbors, and selects a Skyline with rank-based minimum-degree
completion. Queries use weighted-distance graph search with an `ef` search budget.

[Experimental weights](docs/WEIGHTS.md) · [Data formats](docs/FORMAT.md) ·
[Implementation notes](docs/IMPLEMENTATION.md)

This repository contains the core implementation, a synthetic example and a
small validation suite.

## Build

Requirements: Linux, a C++17 compiler with OpenMP (tested with GCC), and Make.
The example generator and tests use Python 3.8+ and its standard library.

```bash
make -j2
make test
```

## Run a complete example

The example has four components: a vector with L2 distance, a scalar with L1
distance, a string with edit distance, and a token set with Jaccard distance.
The generator creates separate database and query objects.

```bash
python3 examples/generate.py --output output/example

HJG_THREADS=12 ./sjg build output/example/database.txt \
  output/example/graph.ind output/example/graph.sjg

EXTERNAL_QUERY_DATA=output/example/query_objects.txt \
  EF=128 OUT_DIST=1 OUT=output/example/answers.txt \
  ./sjg search output/example/database.txt output/example/graph.sjg \
  output/example/queries.txt

EXTERNAL_QUERY_DATA=output/example/query_objects.txt \
  OUT_DIST=1 OUT=output/example/exact.txt \
  ./sjg brute output/example/database.txt output/example/queries.txt
```

Change the first line of `queries.txt` to change the weights; reuse the same graph.
`brute` computes exact neighbors by scanning the database. With
`EXTERNAL_QUERY_DATA`, query IDs refer to the separate query file. Otherwise,
they refer to database objects, including self-matches.

## Main parameters

| Variable | Default | Meaning |
|---|---:|---|
| `HJG_THREADS` | 1 | Build threads |
| `SJG_INDEX_K` | 15 | Per-component target neighbor count |
| `ND_K` | K + 1 | NN-Descent candidate-pool size |
| `ND_EMIT` | K | Neighbors emitted per component |
| `ND_ITERS` | 12 | Maximum NN-Descent iterations |
| `SKY_FLOOR` | 16 | Minimum-degree target, limited by available candidates |
| `EF` | 128 | Search candidate/result budget; raised to k if necessary |
| `EF_LIST` | unset | Comma-separated budgets evaluated in one process |
| `SJG_EDGE_BOUNDS` | 0 | Set to `1` to enable triangle-inequality edge filtering in EF search |
| `OUT` | `ans-tbg.txt` | Answer file; `EF_LIST` appends `.<ef>` |
| `OUT_DIST` | unset | Include individual returned distances |
| `EXTERNAL_QUERY_DATA` | unset | Separate query-object payload file |

Skyline selection is enabled by default. Query execution is single-threaded.
The implementation uses deterministic per-node/per-component initialization;
thread scheduling can still affect parallel build execution. Use one build
thread when checking exact reproducibility.

See [data and output formats](docs/FORMAT.md) and
[implementation notes](docs/IMPLEMENTATION.md).

## Experimental query weights

The [weight reference](docs/WEIGHTS.md) lists the component order and all seven
query-weight profiles for each of the six benchmark datasets. Download the
[full-precision CSV](docs/weights.csv) to reuse the exact experimental values.
The same SJG graph serves all seven profiles.

## Repository layout

```text
src/          construction, metrics, search, I/O and command-line entry
examples/     synthetic example generator
tests/        independent distance checks and end-to-end smoke test
docs/         formats, experimental weights and implementation notes
```

## Search options

EF search evaluates high-weight components first and stops distance evaluation
when the accumulated value exceeds the candidate threshold.
Set `SJG_EDGE_BOUNDS=1` to additionally filter edges using triangle-inequality
bounds. The filter uses component distances stored in the `.sjg` index.

## Baseline source repositories

The experimental baselines use the following upstream implementations:

| Baseline | Upstream repository |
|---|---|
| C-HNSW | [nmslib/hnswlib](https://github.com/nmslib/hnswlib) |
| C-NSG | [ZJULearning/nsg](https://github.com/ZJULearning/nsg) |
| C-Vamana | [cmuparlay/ParlayANN](https://github.com/cmuparlay/ParlayANN) |
| HJG | [ZJU-DAILY/HJG](https://github.com/ZJU-DAILY/HJG) |

The C- prefix denotes our composition of per-component indexes, with candidate
merging and weighted-distance reranking. Our experimental code provides the
composition wrappers and multi-metric adapters. C-Vamana uses the ParlayANN
implementation of Vamana. Baseline source code is available from the repositories
above under their respective licenses.

## License

[MIT](LICENSE).
