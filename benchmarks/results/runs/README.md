# Raw benchmark output

Unedited output of the benchmark runs behind the reported numbers, kept so
each can be checked against what the script printed.
Files are named by date and by the run; most carry their commit, machine and
parameters in their first lines.

| file (after the date prefix) | script | result |
|---|---|---|
| `gap_depth20.txt`, `gap_depth40.txt` | `benchmark_heuristic_gap.py` | at the deployed cap the objective gap stays |
| `trio_truth.txt` | `build_trio_truth.py` | a truth set no phaser touched |
| `oq7_small_dp.txt`, `oq7_small_bf.txt`, `oq7_scale_h*.txt` | `benchmark_genotype_errors.py` | OQ-7 at scale |
| `adaptive_signal.txt` | `benchmark_adaptive_signal.py` | same |
| `convention_cost.txt` | `benchmark_convention_cost.py` | what the span convention costs |
| `scale.txt`, `scale_merge.txt` | `benchmark_selection_at_scale.py` | milestone 2; beyond week 10 |
| `adaptive_scale.txt`, `adaptive_scale_posterior.txt` | `benchmark_adaptive_at_scale.py` | beyond week 10 |
| `corr_none.txt`, `corr_systematic.txt`, `corr_burst.txt`, `corr_both.txt` | `benchmark_selection_at_scale.py` | merge-then-exact under unreported errors |
| `adaptive_scale_stitched.txt`, `adaptive_scale_shuffled.txt`, `adaptive_scale_control.txt`, `adaptive_signal_control.txt` | `benchmark_adaptive_at_scale.py`, `benchmark_adaptive_signal.py` | M4: stitching, and the control |
| `m4_window_*.txt`, `m4_iter_*.txt`, `m4_boost_*.txt`, `m4_merge.txt` | `benchmark_adaptive_at_scale.py` | M4 sensitivity, and merging plus the loop |
| `heldout_window_25.txt`, `heldout_window_50.txt` | `benchmark_adaptive_at_scale.py` | the stitching window, checked on held-out instances |
| `both_original_stitched.txt`, `both_merged_stitched.txt` | `benchmark_adaptive_at_scale.py` | the best pipeline under unreported errors |
| `wh_gap_depth20.txt`, `wh_gap_depth40.txt`, `wh_scale.txt`, `wh_corr_none.txt`, `wh_corr_both.txt` | `benchmark_heuristic_gap.py`, `benchmark_selection_at_scale.py` | the headline was measured on the paper's selector; block length joins the tables |
| `whatshap_score_rule.txt` | `diagnose_whatshap_score_rule.py` | one line of WhatsHap fragments phase blocks |
| `resources.txt` | `benchmark_resources.py` | runtime and peak memory |
| `m4p_*.txt` | `benchmark_adaptive_at_scale.py --merge` | M4: the best pipeline's own sensitivity |
| `pedigree.txt` | `benchmark_pedigree.py` | M4: in pedigree mode the parents make read selection almost irrelevant |
| `hg002_score.txt`, `hg002_score_posthoc_mask.txt`, `hg002_readlists.txt`, `hg002_jobs.txt` | `score_hg002.py`, `compare_whatshap_reads.py`, `hg002_pipelines.py` | M5, weeks 19 to 21: the registered HG002 run (real data) |
| `2026-09-27_oq15_stitched.txt`, `2026-09-27_oq15_oriented.txt` | `benchmark_adaptive_at_scale.py` | OQ-15 repaired |
| `2026-09-27_hg002_chr20_posthoc_new_arms.txt` | `score_hg002.py --exclude 10793-10820` | the new arms on chromosome 20, post hoc (real data) |
| `2026-09-27_trio_budget.txt` | `benchmark_trio_budget.py` | T3: trio budgets |
| `2026-09-27_genome_truth_and_depth.txt` | `build_trio_truth.py`, `genome_run.py` | the genome-wide truth and read depths (real data) |
| `hg002_describe.txt`, `hg002_smoke.txt` | `describe_fragments.py`, the external tools and `extract_fragments.py` | M5, week 18 (real data: HG002 chromosome 20) |
| `het_old_scoring.txt`, `het_spanned.txt`, `confound_spanned.txt` | `benchmark_heterozygosity.py`, `benchmark_margin_confound.py` | corrected: block breaks were scored |

The milestone-1 grid is `../phase_diagram.csv`, which carries its own
provenance. All of it is simulated data except the `hg002_` files, which come from HG002 reads.
