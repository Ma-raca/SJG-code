# Implementation notes

The six C++ source/header files come from the SJG research implementation.
This repository provides the construction and search code, a command-line
interface, a synthetic example and validation tests.

## Construction

SJG generates nearest-neighbor candidates with per-component NN-Descent,
unites forward and reverse candidates, and selects neighbors using
component-distance dominance. Rank-sum completion adds candidates to reach
the minimum-degree target. The score sums each candidate's component-wise ranks.

The command-line defaults enable Skyline selection, a minimum-degree target
of 16, one construction thread and EF 128. The parameter table in the
[README](../README.md#main-parameters) lists the available controls.

## Search

EF search traverses the joint graph using the query's weighted distance.
It evaluates high-weight components first and stops distance evaluation
when the accumulated value exceeds the candidate threshold.

`SJG_EDGE_BOUNDS` defaults to `0`. Set it to `1` to enable triangle-inequality
edge filtering in EF search. The native `.sjg` format stores component distances
on every edge for this filter. Index storage is the same for both settings.

## External queries and formats

After loading and reordering the graph, the external-query adapter appends
query payloads to the distance-access storage. Graph vertices and result IDs
refer to database objects; external queries use the separate payload range.

The implementation supports native `.sjg` indexes, legacy `.ind` indexes and
the earlier EPSI search path. `SJG_EDGE_BOUNDS` controls EF search;
EPSI uses its own bound propagation. Distance and file conventions follow
the research implementation of HJG. Construction uses the NN-Descent
candidate-generation approach.

## Validation

`make test` checks the full pipeline against an independent Python distance
oracle. Cases cover external queries, mixed L1/L2/edit/Jaccard components,
zero and nonuniform weights, and triangle-inequality filtering with both settings.
