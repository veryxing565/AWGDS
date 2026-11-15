# AWGDS

```text
  ###   #   #   ####  ####    ####
 #   #  #   #  #      #   #  #
 #####  # # #  #  ##  #   #   ###
 #   #  ## ##  #   #  #   #      #
 #   #  #   #   ###   ####   ####
```

**Ancient Whole-Genome Duplication Search**  
A command-line workflow for ranking internal branches by retained gene duplication and checking for chromosome-scale evidence of ancient whole-genome duplication (WGD).

**Python ≥ 3.10 · Linux or WSL for the full analysis · One ranked result table**

[Installation](#installation) · [Quick start](#quick-start) · [Outputs](#outputs) · [Input formats](#input-formats) · [Parameters](#parameters) · [Methods](#methods)

AWGDS maps within-species A/B gene pairs to OrthoFinder duplication events, ranks branches with the Abnormal Duplication Index (ADI), and examines the recommended branches for recurring outgroup-anchored chromosome links. Existing OrthoFinder outputs, dated trees, posterior traces and AWGDS branch pairs can be reused.

## Installation

Download and unpack the source distribution or open the repository directory. Use Python 3.10 or later.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
awgds --version
awgds --help
```

For a downloaded release wheel, run `python -m pip install ./awgds-0.3.0-py3-none-any.whl` in your selected Python environment. `python -m awgds` is equivalent to the `awgds` command.

| Analysis | External executables |
| --- | --- |
| Fresh orthogroups and species tree | OrthoFinder, DIAMOND and the helper tools configured for that OrthoFinder installation |
| Pair similarity with reused OrthoFinder results | DIAMOND. A compatible `--diamond-hits` table can supply this evidence |
| Automatic dating with `--divergence-times` | IQ-TREE 3 with `--dating mcmctree`, PAML MCMCtree, and MAFFT for automatic alignments |
| Revalidation from saved pairs | Python dependencies installed with AWGDS |

Put the relevant executables on your `PATH` or pass their locations as [parameters](#parameters). Start OrthoFinder from the environment containing its helper programs. [OrthoFinder](https://orthofinder.github.io/OrthoFinder/tutorials/beginner-tutorial/) · [DIAMOND](https://github.com/bbuchfink/diamond) · [IQ-TREE dating](https://iqtree.github.io/doc/Dating) · [PAML](https://github.com/abacus-gene/paml) · [MAFFT](https://mafft.cbrc.jp/alignment/software/)

## Quick start

### 1. Run the whole analysis

```bash
awgds genomes/ annotations/ -o results/ --threads 16 \
  --chromosome-info chromosomes.txt
```

Use FASTA genomes and matching GFF3 files with CDS features. Each output directory must be new or empty. With a chromosome-only FASTA and GFF3 collection, omit `--chromosome-info`. With no dating input, AWGDS uses the original OrthoFinder branch lengths. Open `results/summary.tsv` when the run completes.

### 2. Reuse OrthoFinder and a dated tree

```bash
awgds genomes/ annotations/ -o dated_results/ --threads 16 \
  --orthofinder-results Results_Nov09/ \
  --chromosome-info chromosomes.txt \
  --time-tree clock=chronogram.nwk --time-unit Myr
```

For two clocks, repeat `--time-tree` and supply a common `--time-unit` or one unit per input. To date the species tree from calibration points or intervals, use `--divergence-times divergence_times.txt`. See [Dating inputs](#dating-inputs).

### 3. Revalidate the recommended branches

```bash
awgds validate results/ -o revalidated/ --threads 16 --wgd-local always
```

Point `validate` to the **full run directory**. It reads the saved ADI scores, labelled species tree and branch pairs, and uses the recorded genome, GFF3 and OrthoFinder membership inputs. It produces the same unified `summary.tsv` with the original ADI ranking and the current WGD evidence. AWGDS 0.1 and 0.2 full runs are supported. You can supply replacement paths when inputs from the same genome release have moved.

## Outputs

### The main result

**`summary.tsv`** contains one row per scored internal branch in ADI rank order. It includes branches below the recommendation threshold, with their WGD evaluation status.

| Columns | Contents |
| --- | --- |
| `rank`, `node`, `parent` | Order and labels from the rooted species tree |
| `density_unit`, `clock_density_<name>` | One median duplicate-family density column per clock |
| `adi_raw`, `adi_percent`, `recommended` | Raw ADI, its relative 0–100 scale and the recommendation flag |
| `global_ratio`, `neighbor_ratio` | Conservative ratios against the whole-tree and neighboring-branch backgrounds |
| `wgd_status`, `wgd_selected_references` | Validation outcome and the selected outgroup species |
| `wgd_tested_target_reference_combinations`, `wgd_positive_target_reference_combinations`, `wgd_positive_pairs` | Test counts and the exact positive combinations, formatted as `target\|reference` entries separated by semicolons |
| `wgd_reference_coverage` | Shared versus testable chromosomes for each reference, for example `Taxon_C=3/3;Taxon_E=2/4` |
| `wgd_best_reference`, `wgd_shared_reference_chromosomes`, `wgd_testable_reference_chromosomes`, `wgd_shared_chromosome_fraction` | The strongest reference and its chromosome evidence |
| `wgd_two_sided_references`, `wgd_left_eligible_targets`, `wgd_right_eligible_targets` | Availability of comparisons across the daughter split |
| `wgd_local_check`, `wgd_local_matches`, `wgd_local_target_pairs` | Local-chain status, match count and matched daughter-side species pairs |

`shared_macrosynteny` means the cross-daughter chromosome evidence meets the configured thresholds. `no_broad_linkage_detected` means an eligible comparison fell below the coverage criterion. `insufficient_data` records limited chromosome or anchor coverage. `not_recommended`, `not_scored` and `validation_disabled` identify branches outside the validation run. `pending` marks a saved ADI checkpoint while validation is in progress or incomplete. A local check with zero matches has `0` in its match column. A skipped local check has `NA`.

### Retained files

```text
results/
|-- summary.tsv                         # Unified final result
|-- run.json                            # Parameters, provenance, thresholds and status
|-- logs/
|   |-- run.log                         # Banner, numbered stages and output paths
|   `-- tools/                          # Full logs for external programs, when run
`-- work/
    |-- inputs/                         # Genome/GFF paths, extraction report and input maps
    |-- orthofinder/                     # Complete OrthoFinder results from a fresh run
    |-- trees/                          # Rooted tree, archived time trees and dating runs
    |-- similarity/                     # Reusable DIAMOND hits and optional raw search files
    |-- branches/
    |   |-- pairs/                      # One compressed gene-pair table per node
    |   |-- branch_taxa.tsv              # Taxa on the two sides of each split
    |   |-- branch_counts.tsv           # Counts by filtering step
    |   `-- events.tsv.gz               # Event-level copy and pair evidence
    |-- adi/                            # Family panel, observations and per-clock intervals
    `-- evidence/                       # Detailed WGD tests for recommended branches
```

| Path | What to keep it for |
| --- | --- |
| `run.json` | Run status, versions, options, external commands, tree sources and the ADI recommendation threshold. Check for `"status": "completed"` before using a final run |
| `work/inputs/input_files.tsv`, `extraction.tsv` and optional maps | Genome and annotation provenance. Fresh OrthoFinder inputs are available in `work/inputs/proteins/` |
| `work/orthofinder/` | OrthoFinder's native result directory, including the original labels, memberships and duplication calls. In reuse mode the original result path is recorded in `run.json` |
| `work/trees/species_tree.nwk` | Labelled rooted species tree used for branch assignments |
| `work/trees/imported/`, `work/trees/mcmctree/` | Saved time-tree inputs, compressed imported posterior traces and diagnostics, or generated dated trees, calibration settings, alignments and chains |
| `work/similarity/diamond_hits.tsv.gz` | Candidate hits suitable for `--diamond-hits` with the same eligible genes and search settings. `--keep-search` keeps candidate FASTAs, identifier maps, databases and native hits under `work/similarity/raw/` |
| `work/branches/pairs/<node>.pairs.tsv.gz` | A/B IDs, species, orthogroup, gene-tree split, GFF3 coordinates, replicons and hit metrics. `used_in_adi=1` marks pairs from contributing events in the common ADI panel |
| `work/adi/families.tsv`, `family_species_pairs.tsv.gz`, `clock_densities.tsv` | Common family panel, deduplicated ADI observations, branch-length and density intervals, and comparison backgrounds |
| `work/evidence/selection.tsv`, `chromosomes.tsv` | Reasons for selecting targets and outgroups, and the sequence-to-chromosome assignments used |
| `work/evidence/target_reference.tsv`, `chromosome_links.tsv.gz` | Comparison-level anchor counts and chromosome-pair tests with max-adjusted permutation P-values |
| `work/evidence/reference_summary.tsv`, `shared_linkage.tsv` | Eligible chromosomes per reference and matched A/B events across daughter lineages |
| `work/evidence/local_chains.tsv` | Local chain matches across branches, with node, species, sequence ranges and event IDs. Created when the local check runs |
| `work/evidence/matched_families.tsv` | Optional family-count-matched comparisons across branches |

Revalidation writes a new unified summary and evidence tables. Its `run.json` links to the original full run for the branch pairs and tree files.
The evidence directory is created when at least one branch is evaluated. Dating and raw search directories appear when their corresponding options are used.

## Input formats

### Genome FASTA and GFF3

AWGDS searches the two input directories recursively. The case-sensitive species key is the filename prefix before the first dot. For example, `Taxon_A.softmasked.fna.gz` matches `Taxon_A.longest.gff3`. Species keys start with a letter or digit and may contain letters, digits, underscores and hyphens. Provide one genome and one GFF3 per key, with at least four species for the branch comparison.

FASTA extensions are `.fna`, `.fa`, `.fasta` and `.fas`. Annotation extensions are `.gff3` and `.gff`. Both accept a final `.gz`. Use the FASTA header's first token as the GFF3 sequence ID. The GFF3 should contain phased CDS records and gene/transcript relationships. Here is one record of each format.

```text
>chr01 Taxon_A
ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTTAA
```

```text
##gff-version 3
chr01	Example	gene	1	66	.	+	.	ID=g1
chr01	Example	mRNA	1	66	.	+	.	ID=g1.t1;Parent=g1
chr01	Example	CDS	1	66	.	+	0	Parent=g1.t1
```

For a fresh run, the longest usable CDS translation per gene is selected. Reuse assemblies and annotations from the same releases used to generate the OrthoFinder proteins. [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

### Chromosome assignments

`--chromosome-info` accepts a whitespace-delimited whitelist with an optional header. Enter a species key or genome basename and the exact FASTA sequence ID. An optional third column gives the chromosome label and can join components assigned to the same chromosome.

```text
species	seqid	chromosome
Taxon_A	CM012345.1	chr01
Taxon_A	CM012346.1	chr02
Taxon_B	chr01	chr01
Taxon_B	chr02	chr02
```

In whitelist mode, omitted sequences remain unassigned. Without a whitelist or replicon map, AWGDS assumes all supplied FASTA sequences are chromosomes. [Chromosome whitelist example](examples/chromosomes.txt) · [Legacy replicon map example](examples/replicon_map.tsv)

### Reusable orthogroups and protein hits

`--orthofinder-results` points to an OrthoFinder result directory with the labelled species tree, duplication calls and orthogroup tables. These paths are standard.

```text
Results_Nov09/
|-- Species_Tree/SpeciesTree_rooted_node_labels.txt
|-- Gene_Duplication_Events/Duplications.tsv
|-- Orthogroups/Orthogroups.tsv
|-- Orthogroups/Orthogroups.GeneCount.tsv
`-- WorkingDirectory/                   # Archived protein FASTAs, when available
```

The species-tree tips and orthogroup columns must match the genome keys. Chromosome verification uses orthogroup membership to anchor outgroup singleton genes. AWGDS can reuse the archived protein sequences in a complete `WorkingDirectory`. Pass `--gene-id-map` to resolve individual OrthoFinder IDs against your annotation.

```text
species	orthofinder_id	gff_id
Taxon_A	old_transcript_100	g1.t1
```

`--diamond-hits` accepts AWGDS's normalized hit export. `species`, `gene1`, `gene2`, `evalue` and `bitscore` are required. Identity and query/subject coverage are percentages and support the corresponding filters.

```text
species	gene1	gene2	evalue	bitscore	pident	qcov	scov
Taxon_A	g1.t1	g2.t1	1e-30	142	61.5	93	91
```

[ID map example](examples/gene_id_map.tsv)

### Dating inputs

**Fixed trees.** `--time-tree` accepts a dated, rooted Newick tree, a FigTree/NEXUS tree or the dated mean tree in an MCMCtree report. Supply `--time-unit` as `years`, `kyr`, `Myr`, `100Myr` or `Gyr`. Target species must be present. After extra tips are pruned, the rooted clades must agree with the OrthoFinder tree and the projected time tree must be ultrametric.

```newick
((Taxon_A:5,Taxon_B:5):10,(Taxon_C:5,Taxon_D:5):10);
```

**Posterior draws.** Provide the matching MCMCtree report and its comma-separated traces to propagate dating uncertainty. The report must include the numbered tree and the traces must contain its age columns, such as `t_n9`. Supply the report's original length unit.

```bash
awgds genomes/ annotations/ -o posterior_run/ \
  --orthofinder-results Results_Nov09/ \
  --mcmc-posterior independent dating/clock2_chain1/mcmctree.out \
    dating/clock2_chain1/mcmc.txt,dating/clock2_chain2/mcmc.txt \
  --time-unit 100Myr
```

When mixing fixed trees and posterior bundles, list repeatable units in the order of the fixed trees followed by the posterior bundles. A single unit applies to all inputs.

**Calibrations.** `--divergence-times` takes a tab-separated table. In each row, `species_a` and `species_b` define the calibrated MRCA. Supply both `min_mya` and `max_mya` or a positive `age_mya`. By default, a point age becomes a soft interval of ±20%. The root needs either a calibration row or `--root-max-mya`. An optional `source` column records where an age came from.

```text
species_a	species_b	age_mya	min_mya	max_mya	source
Taxon_A	Taxon_D		170	260	Illustrative root range
Taxon_A	Taxon_B	100			Illustrative point age
```

AWGDS uses complete single-copy groups for the MAFFT and IQ-TREE alignment workflow by default. `--dating-alignment` accepts a pre-aligned protein FASTA or relaxed PHYLIP with one record per sampled species. Review the saved diagnostics for MCMCtree chains before interpreting inferred durations. [Calibration example](examples/divergence_times.txt)

## Parameters

Run `awgds --help` for the full workflow or `awgds validate --help` for validation of a saved run.

<details>
<summary><strong>Inputs and execution</strong></summary>

| Argument | Default | Purpose |
| --- | --- | --- |
| `GENOMES`, `GFFS` | Required | Folders containing matched genome and annotation files |
| `RUN` | Required for `validate` | Saved full-run directory |
| `-o`, `--outdir DIR` | Required | New or empty output directory |
| `-t`, `--threads N` | Up to 8 | Threads for each parallel stage |
| `-h`, `--help` | — | Display command help |
| `--version` | — | Display the AWGDS version |
| `--genomes DIR`, `--annotations DIR` | Saved source paths | In `validate`, supply both to repoint to matching genome and GFF3 folders |
| `--orthofinder-results DIR` | Fresh OrthoFinder run | Reuse complete results, or override the saved location in `validate` |
| `--species-tree TREE` | Inferred by OrthoFinder | Supply a rooted, bifurcating species tree for a fresh run |
| `--orthofinder-arg TOKEN` | None | Repeat to pass tokens to fresh OrthoFinder. Example `--orthofinder-arg=-ST --orthofinder-arg=astral-pro` |
| `--orthofinder EXE`, `--diamond EXE` | `orthofinder`, `diamond` | Select executable names or paths |
| `--protein-id {transcript,protein,gene}` | `transcript` | IDs in fresh protein FASTAs |
| `--prefix-protein-ids` | Off | Prefix fresh FASTA protein IDs with the species key |
| `--genetic-code INT` | `1` | Default NCBI translation table. CDS-level `transl_table` takes precedence |
| `--min-protein-aa N` | `20` | Minimum usable protein length |
| `--retranslate` | Off | Derive proteins from CDS when reusing OrthoFinder results |
| `--gene-id-map TSV` | None | Resolve alternate gene IDs, or override the recorded map in `validate` |
| `--allow-unmapped-genes` | Off | Permit saved IDs with missing coordinates when enabled filters allow them |

</details>

<details>
<summary><strong>Gene-pair filters and chromosomes</strong></summary>

| Argument | Default | Purpose |
| --- | --- | --- |
| `--same-og` / `--no-same-og` | On | Require membership in the orthogroup named by the duplication record |
| `--cross-replicon` / `--no-cross-replicon` | On | Require the two genes to lie on different assigned chromosome or sequence IDs during candidate filtering |
| `--require-diamond` / `--no-require-diamond` | On | Require the branch-pair similarity hit. Fresh OrthoFinder still uses DIAMOND for its own search |
| `--cross-daughters` / `--no-cross-daughters` | On | Export only events with accepted A/B pairs on both daughter sides. One-sided events remain outside ADI |
| `--min-dup-support FRAC` | `0` | Minimum OrthoFinder duplication support on the 0–1 scale |
| `--family-occupancy FRAC` | `1` | Minimum fraction of sampled species represented in the common ADI group panel |
| `--chromosome-info TXT` | None | Use an explicit chromosome whitelist and optional chromosome labels |
| `--replicon-map TSV` | None | Use the legacy partial map with `species`, `seqid` and `replicon` columns |
| `--assume-chromosomes` | Off | In `validate`, use all input sequence IDs in place of the saved chromosome map |
| `--diamond-hits TSV` | New DIAMOND search | Read saved candidate hits |
| `--evalue E` | `1e-5` | Maximum DIAMOND E-value |
| `--min-identity PCT` | `0` | Minimum DIAMOND identity percent |
| `--min-coverage PCT` | `0` | Minimum percent coverage of each protein |
| `--max-targets N` | `0` | Cap hits per query. Zero retains all passing targets |
| `--rbh` | Off | Keep reciprocal best eligible pairs, ranked by bit score, E-value, then ID |
| `--keep-search` | Off | Retain native DIAMOND inputs, ID maps, databases and results |

An updated chromosome whitelist in `validate` operates on the candidate pairs already present in the saved run. During a full run, chromosome assignments apply to both candidate filtering and validation. The macro test always uses pairs on different chromosomes. Same-chromosome pairs can enter the fallback local check when the original full run includes them via `--no-cross-replicon`.

</details>

<details>
<summary><strong>Time trees and ADI</strong></summary>

| Argument | Default | Purpose |
| --- | --- | --- |
| `--time-tree [NAME=]TREE` | None | Dated Newick, NEXUS/FigTree or MCMCtree mean tree. Repeatable |
| `--mcmc-posterior NAME REPORT TRACE[,TRACE]` | None | Numbered MCMCtree report with one or more traces. Repeatable |
| `--time-unit UNIT` | Required for imported dates | `years`, `kyr`, `Myr`, `100Myr` or `Gyr` |
| `--divergence-times TXT` | None | Calibrate the species tree from MRCA ages |
| `--point-age-width FRAC` | `0.2` | Fraction of a point age on each side of its soft interval |
| `--root-max-mya AGE` | None | Root age upper bound when the table has no root calibration |
| `--dating-alignment FILE` | Auto-align single-copy families | Use a supplied protein alignment |
| `--dating-max-ogs N` | `0` | Select at most this many single-copy orthogroups. Zero uses all |
| `--dating-model MODEL` | `LG+F+G4` | IQ-TREE protein model for dating |
| `--dating-clock {both,independent,autocorrelated}` | `both` | Clock models for calibration, with two chains per clock |
| `--mcmc-burnin N` | `100000` | Burn-in per chain |
| `--mcmc-samples N` | `18000` | Saved samples per chain |
| `--mcmc-sample-every N` | `20` | Iterations between saved MCMC samples |
| `--iqtree EXE`, `--mafft EXE`, `--mcmctree EXE` | `iqtree3`, `mafft`, `mcmctree` | Executable names or paths for dating |
| `--bootstrap N` | `12000` | ADI family replicates |
| `--seed N` | `410611` | Reproducible random seed. Revalidation inherits the source seed by default |

</details>

<details>
<summary><strong>WGD validation</strong></summary>

| Argument | Default | Purpose |
| --- | --- | --- |
| `--wgd-validation` / `--no-wgd-validation` | On | Run or skip chromosome validation after ranking |
| `--wgd-permutations N` | `1999` | Permutations per target/reference comparison |
| `--wgd-alpha P` | `0.05` | Strict threshold for within-comparison max-adjusted P |
| `--wgd-min-anchor-families N` | `20` | Minimum shared families for a target/reference test |
| `--wgd-min-link-families N` | `5` | Minimum sum of family weights for a chromosome link |
| `--wgd-min-shared-events N` | `3` | Minimum common events with coherent A/B orientation across daughter sides |
| `--wgd-min-shared-chromosomes N` | `2` | Minimum matching reference chromosomes for shared macrosynteny |
| `--wgd-min-placed-fraction FRAC` | `0.5` | Minimum fraction of coding loci on assigned chromosomes |
| `--wgd-max-references N` | `4` | Maximum selected outgroups per branch. Zero uses all eligible outgroups |
| `--wgd-max-outgroup-depth N` | `0` | Limit successive sister clades searched. Zero searches all depths |
| `--wgd-exclude-outgroup SPECIES` | None | Exclude a specified reference species. Repeatable |
| `--wgd-local {fallback,always,off}` | `fallback` | Run local chains after insufficient macro evidence, on every recommended branch, or skip them |
| `--wgd-local-gap N` | `25` | Maximum annotated coding loci between adjacent ordered anchors |
| `--wgd-local-min-chain N` | `5` | Minimum distinct families per local chain |
| `--wgd-matched-draws N` | `0` | Enable family-count-matched resampling for the same comparison on multiple branches |
| `--wgd-matched-permutations N` | `499` | Permutations per matched resample |

</details>

## Methods

1. **Identify branch pairs.** AWGDS reads the A/B clades of each OrthoFinder duplication record and finds pairs within the same sampled descendant species. By default, pairs must satisfy orthogroup, DIAMOND and chromosome filters and have support from both daughter sides at the event level.
2. **Rank retained duplication.** Each family contributes at most once for a particular cross-daughter species pair. Contributions are averaged over possible species pairs for the branch and divided by branch length, per 1,000 common families. AWGDS bootstraps families and optionally samples joint posterior branch durations. For each reference background, ADI uses the minimum across clocks of the branch's lower density limit divided by the comparison background. The larger of the global and neighbor ratios becomes the score. Candidates are recommended above `Q3 + 1.5 × IQR` on raw ADI.
3. **Identify chromosome evidence.** Outgroups are drawn from successive sister clades outside the focal branch, prioritizing references with testable descendants on both sides. A singleton reference gene anchors its family and each family contributes a total weight of one across possible target chromosome pairs. A max-statistic permutation tests chromosome associations within each target/reference comparison. Shared evidence requires common duplication events and coherent A/B orientation on the two daughter sides. A local-chain check adds evidence for ordered segments when useful.

ADI describes enrichment among retained, detectable families. The permutation P-values quantify chromosome association under the selected comparison, and the branch-level evidence status describes agreement across daughter lineages. Copy loss, chromosome rearrangement, family reconstruction, assembly and species sampling all affect detection.

## Troubleshooting

| Message or result | What to check |
| --- | --- |
| Species names differ | FASTA and GFF3 filename prefixes, tree tip labels and orthogroup column names, including case |
| Missing or non-positive branch lengths | Use a compatible rooted tree with positive internal branch lengths or supply explicit time inputs |
| No recommendable ADI scores | Check the shared family panel, candidate filters and retained A/B pairs in `work/branches/branch_counts.tsv` |
| Limited chromosome evidence | Check the whitelist and `work/evidence/selection.tsv` for eligible outgroups and both daughter sides |
| An external command failed | Open the command-specific log under `logs/tools/` and check the executable and its dependencies |

To run the regression tests from the source tree after installation:

```bash
python -m unittest discover -s tests -v
```

## Citation

Record the AWGDS version or commit, the command line and the saved `run.json` when reporting results. Cite the methods used for the underlying analyses. [OrthoFinder](https://doi.org/10.1186/s13059-019-1832-y) · [DIAMOND](https://doi.org/10.1038/s41592-021-01101-x) · [Comparative macrosynteny](https://pmc.ncbi.nlm.nih.gov/articles/PMC7145513/) · [NIST boxplot fences](https://www.itl.nist.gov/div898/handbook/prc/section1/prc16.htm)
