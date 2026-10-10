# shinyNTApro

shinyNTApro is an R/Shiny application for analysing NanoSight nanoparticle tracking analysis (NTA) data. It combines data import, interactive plots, statistical comparisons, quality control and downloadable HTML reports.

## Features

- Import supported NanoSight NTA 3.4 summary CSV files.
- Visualise particle concentration, size distributions and sample comparisons.
- Compare independent groups using parametric or nonparametric tests, with post-hoc comparisons and effect sizes.
- Monitor measurement performance using Levey–Jennings plots and precision summaries.
- Export plots, statistical results and HTML reports.

## Getting started

Use the [online application](https://andrewlai.shinyapps.io/shinyNTApro/), or run locally:

1. Download this repository.
2. Install R and the packages listed at the beginning of `app.R`.
3. Open R in the repository folder and run:

```r
shiny::runApp()
```

Upload a supported summary CSV, assign sample groups and select the analysis you need. Use a local installation for persistent QC logging.

Automatic test selection provides guidance; users remain responsible for choosing analyses appropriate to their study design.

## Feedback

Report bugs or suggest improvements through [GitHub Issues](https://github.com/andrewlai/shinyNTApro/issues).
