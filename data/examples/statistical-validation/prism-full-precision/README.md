# Prism validation imports and review

Reviewed `Prism_check.prism`, saved by GraphPad Prism 11.1.0.224 on 24 September 2026. The original project was not changed.

## Findings from the saved project

The project contains two input datasets and two analysis configurations. Only the ordinary ANOVA has populated saved numerical results. The Welch t-test is marked `isCalculated: false`, with no numerical results in its saved results JSON. Its settings specify an unpaired, two-sided Welch t-test, alpha 0.05 and the B − A direction. Recalculate and save that analysis before using this project as its validation record.

Both input tables contain rounded values. All 60 values in the two-group dataset and all 90 values in the three-group dataset differ from their corresponding source CSV values. The largest absolute difference is 0.00004637304853361 nm. This has little effect on the displayed results, but identical inputs are needed for a precise software comparison.

| ANOVA quantity | Original source data, calculated in R | Saved Prism result, rounded input |
|---|---:|---:|
| F statistic | 272.26536 | 272.269218017937 |
| Numerator / denominator df | 2 / 87 | 2 / 87 |
| P value | approximately 3.564 × 10⁻³⁸ | 3.562023189192 × 10⁻³⁸ |
| Eta squared / R squared | approximately 0.8622395 | 0.862241163743 |

R reproduces the Prism F statistic when given the rounded data extracted from the project. The discrepancy therefore arises from input rounding.

## Tukey comparisons

Prism reports first group minus second group. The shinyNTApro output supplied for this case reports the opposite direction. Keep the group identities unchanged. To express a confidence interval in the opposite direction, negate the estimate and transform `[lower, upper]` to `[-upper, -lower]`.

| Comparison | Saved Prism difference (nm) | Saved Prism simultaneous 95% CI (nm) | Stored adjusted P value |
|---|---:|---|---:|
| A − B | −3.000003666667 | −3.613103016090 to −2.386904317243 | 3.089979383475 × 10⁻¹⁰ |
| A − C | −6.000003666667 | −6.613103016090 to −5.386904317243 | 3.089899447417 × 10⁻¹⁰ |
| B − C | −3.000000000000 | −3.613099349423 to −2.386900650577 | 3.089979383475 × 10⁻¹⁰ |

The saved numerical probabilities are available even though the displayed Prism table uses P < 0.0001. Prism and R both return values near 3.09 × 10⁻¹⁰ here. This agreement does not establish extreme-tail accuracy: for these balanced three-group comparisons, Bonferroni upper bounds on the true Tukey tail probabilities are approximately 5.23 × 10⁻¹⁹ for adjacent groups and 9.80 × 10⁻³⁹ for A versus C. The much larger reported values indicate a numerical tail limitation. Preserve the raw outputs, describe the limitation, and use an appropriately qualified threshold in the supplementary presentation instead of claiming precise validation of these extreme probabilities.

## Importing the full-precision data

Each CSV below has one column per group (A, B and, where applicable, C). Values retain their original decimal strings. Blank cells retain missing observations or pad shorter groups; they must not be replaced by zero. Measurements are mean average diameter in nm.

Import each CSV into a separate Prism Column table, using its first row as column titles. Use the scenario name as the table title. Import the CSV directly rather than copying rounded values from a displayed table. Keep the independent-samples design and select the matching test and multiple-comparison method. Recalculate and save the project after replacing the two existing datasets.

| Scenario import | Expected shinyNTApro Auto route | Post-hoc method |
|---|---|---|
| [normal_equal_3x30.csv](normal_equal_3x30.csv) | Ordinary one-way ANOVA | Tukey HSD |
| [normal_unequal_3x30.csv](normal_unequal_3x30.csv) | Welch ANOVA | Games–Howell |
| [skewed_3x30.csv](skewed_3x30.csv) | Kruskal–Wallis | Dunn with Holm adjustment |
| [normal_2x30.csv](normal_2x30.csv) | Welch t-test | None |
| [skewed_2x30.csv](skewed_2x30.csv) | Wilcoxon rank-sum | None |
| [small_3x3.csv](small_3x3.csv) | Ordinary one-way ANOVA | Tukey HSD |
| [missing_3x3.csv](missing_3x3.csv) | Auto stops: only two usable observations in A | Not applicable |
| [tied_skewed_3x30.csv](tied_skewed_3x30.csv) | Kruskal–Wallis | Dunn with Holm adjustment |
| [unbalanced_equal_10_20_30.csv](unbalanced_equal_10_20_30.csv) | Ordinary one-way ANOVA | Tukey HSD |
| [unbalanced_smallest_widest_10_20_30.csv](unbalanced_smallest_widest_10_20_30.csv) | Welch ANOVA | Games–Howell |
| [unbalanced_largest_widest_10_20_30.csv](unbalanced_largest_widest_10_20_30.csv) | Welch ANOVA | Games–Howell |
| [unbalanced_skewed_10_20_30.csv](unbalanced_skewed_10_20_30.csv) | Kruskal–Wallis | Dunn with Holm adjustment |
| [unbalanced_normal_10_30.csv](unbalanced_normal_10_30.csv) | Welch t-test | None |
| [unbalanced_skewed_10_30.csv](unbalanced_skewed_10_30.csv) | Wilcoxon rank-sum | None |

These are 14 input scenarios, not 14 completed Prism validations. The missing-data scenario tests the application's stopping rule rather than a Prism inferential result.

For the remaining scenarios, match the numerical method as well as the test name. Prism's standard Dunn adjustment differs from Holm: compare unadjusted Dunn results and apply Holm separately if needed. For rank tests, record whether exact or approximate P values, tie handling and continuity correction were used. Do not treat a difference caused by these settings as a software error.

Source files: the fourteen sibling scenario folders, each containing `<scenario>_values.csv`, and `../expected-results.csv`. The imported rows retain within-group source order. A shared row number does not specify paired observations.

References: [GraphPad Dunn calculation](https://www.graphpad.com/guides/prism/latest/statistics/stat_how_the_dunn_method_for_nonpar.htm), [GraphPad Mann–Whitney calculation](https://www.graphpad.com/guides/prism/latest/statistics/how_the_mann-whitney_test_works.htm), [R studentized-range distribution](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Tukey.html).
