# Implementation and provenance

The six C++ source/header files are extracted from the maintained SJG research
implementation. Construction and distance kernels are retained. The EF search
retains the research traversal, with triangle-inequality edge filtering made
opt-in through `SJG_EDGE_BOUNDS=1` (disabled by default).

The release adds an external-query adapter already used by the independent-query
experiments. Query payloads are appended to the distance-access storage after
loading/reordering the graph. Database cardinality and graph membership stay
unchanged, so external queries cannot become returned neighbors.

Release-specific CLI defaults enable Skyline selection, minimum degree 16,
one construction thread and EF 128. English CLI help, Make targets, format
documentation, a synthetic generator and tests are added. Internal environment
variable names are preserved to keep experiment commands recognizable.

Construction uses per-component NN-Descent, the union of forward/reverse
candidates, component-distance dominance, and rank-sum minimum-degree completion.
The rank score uses component-wise ranks rather than a fixed weighted distance.

The native index stores per-edge component distances. EF search uses these
distances for triangle-inequality filtering only when `SJG_EDGE_BOUNDS=1`.
The default skips this filtering. Weight-ordered component evaluation and
early termination of accumulated distances remain enabled in either mode.
The switch changes query behavior, not the serialized index layout or size.

Some compatibility utilities remain in the same core files (legacy index import
and the older EPSI search path, whose bound propagation is unchanged by the
EF-only switch). No baseline implementation, third-party source
tree, private dataset, generated index, experiment result or paper source is
included. Construction follows the NN-Descent candidate-generation approach;
legacy distance and file conventions were maintained for comparisons with the
research implementation of HJG. Algorithm names in comments document provenance,
not additional bundled baseline libraries.

`make test` exercises the full pipeline with external queries, mixed L1/L2/edit/
Jaccard components, zero component weights, nonuniform weights and an independent
Python distance oracle. It is a correctness smoke test, not a performance result.
