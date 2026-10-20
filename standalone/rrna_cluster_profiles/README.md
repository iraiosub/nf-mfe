# Multi-species clustered rRNA profiles

Run from any directory with Python, NumPy and Matplotlib:

```bash
python standalone/rrna_cluster_profiles/plot_profiles.py
```

Defaults to the named TSV in `~/Downloads` and the structural assets in
`~/github/pdb-rna-benchmarking`. Override with `--input`, `--reference-repo`,
`--out`, `--bin-size` or `--arc-min`. Use `--weight rows` for equally weighted
cluster rows instead of the default `count` support weights.

Produces 18S/28S PNG and PDF panels for native human, lifted mouse, lifted rat,
merged sources, and a strict merged subset requiring both arms to have human
coordinates. The source comparison PDF has a page each for 18S and 28S, with human,
lifted mouse, lifted rat, merged, and strict merged shown as separate row facets.
Expansion segments repeat beside every facet; each coverage facet uses its own x
range. Per-RNA PNG/PDF versions are also written. Separate full
source panels use their own scales, with independent coverage ranges, square structural/contact maps and
vertical profiles aligned to the map height.

Assumptions explicitly requested for this table: human `rDNA` and `rRNA_45S`
share the common human U13369 coordinate system; all source and target intervals
are 1-based inclusive. Human arms use original genomic coordinates. Rodent arms
use successful `LO_target_*` coordinates exclusively. Failed mappings are never
replaced by rodent coordinates. 18S/28S boundaries and physical distance maps
come from the existing 4V6X reference assets, with low model coverage masked.
Only complete arms inside a mature RNA enter that RNA's profiles. All rDNA/rRNA_45S partners outside the within-RNA category are grouped
as Other rRNA, including positions outside 18S/28S and boundary-crossing partners.

Main lifted tracks require the rRNA arm to map; the partner may remain unmapped.
Partner classification uses destination annotation where available, otherwise
source annotation. `merged_both_arms_available` excludes unmapped rodent partners.
Species, proteins and samples are pooled within each displayed source; the input
row is the counting unit and no cross-species cluster deduplication is attempted.

Coverage weights each base by `count`, merging overlapping arms in within-RNA
rows so the row contributes once at each base. Plot values are mean coverage per
bin, without RPM normalization because this is an rRNA-only table. Contact maps
use paired cluster interval midpoints, weighted by the same choice; these are
cluster summaries, not reconstructed read-level contacts. Arcs use off-diagonal
bin pairs with at least `--arc-min` raw weight. Empty data remain empty.

Clippy direct coverage marks rRNA arms carrying a supplied peak annotation;
Clippy partner coverage marks the rRNA footprints whose opposite arm carries a
peak annotation. Neither is an independent peak-density or cDNA-score track.
The direct rRNA track is empty for this input. Peak coordinates cannot safely be
placed on rDNA from non-rRNA partner peak IDs. Original unique peak IDs, regions,
and cDNA scores are exported in `clippy_peaks.tsv`, deduplicated within species.
Blank Clippy fields mean no positive annotation supplied, not proven absence.

`profiles.tsv`, `partner_regions.tsv`, `summary.tsv`, per-source/RNA `.npz`
arrays, and `audit.json` provide underlying numbers, assumptions, mapping counts
and the input SHA-256. Contact arrays store each pair once in the upper triangle;
the displayed map mirrors that triangle. Internal checks require exact merged
coverage/contact sums and agreement between contact totals and within-RNA weights.

`species_aligned.pdf` and `.png` add a compact two-row comparison (18S and
28S): each coverage category contains adjacent narrow Mouse / Rat / Human
tracks aligned to mature RNA positions, with independent coverage ranges and
expansion segments alongside. Existing full source facets remain available.

`interaction_types.species.pdf` and `.png` show partner annotation bar charts
with Human / Mouse lifted / Rat lifted columns and 18S / 28S rows. Counts and
within-facet percentages are labelled, using independent count ranges.

To generate the equally weighted cluster-row version alongside the default:

```bash
python standalone/rrna_cluster_profiles/plot_profiles.py --weight rows --out standalone/rrna_cluster_profiles/results_cluster_rows
```

In this version, each input cluster row contributes one, including to bar-chart
counts and percentages. It does not deduplicate cluster labels across distinct
RNA pairs or species. The support-weighted results remain in `results/`.
