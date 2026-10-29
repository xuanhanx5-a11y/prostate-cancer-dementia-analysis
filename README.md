# Systemic vulnerability and dementia in prostate cancer

Configurable R workflows for cohort description, survival analysis, competing risks, missing-data analysis and latent class analysis.

## Installation

```sh
Rscript environment/install.R
```

Package versions are listed in `environment/tested_package_versions.csv`.

## Configuration and inputs

Place authorised input files in `data/` and configure the relevant analysis:

- Clinical Cox workflow: `code/clinical/cox_config.R`; input `clinical_survival.csv`.
- General survival workflow: `code/general/config.R`; input `analysis.csv`.
- Latent class analysis: `lca_clinical.csv` and `lca_ukb.csv`.

Survival inputs require a unique `id`, positive follow-up `time` and `status`: 0 = censored, 1 = dementia, 2 = death before dementia. Specify the time origin, unit, endpoint definition, exposures, adjustment variables and factor reference levels in the configuration. Time ends at dementia, preceding death or censoring, whichever occurs first.

See [Input and methods](docs/input_and_methods.md) for variable requirements and calculation details.

## Run

From the repository root:

```sh
Rscript code/main.R check
Rscript code/main.R clinical
Rscript code/general/pipeline.R
Rscript code/lca/run_lca.R clinical
Rscript code/lca/run_lca.R ukb
```

`check` verifies syntax and input availability. The clinical command runs baseline summaries, Cox screening, adjusted models, configured subgroup analyses, proportional-hazards diagnostics and an HR forest plot.

Individual clinical stages are `recode`, `baseline`, `screening`, `multivariable`, `subgroups` and `screen_plot`. Run `screening` before `screen_plot`. The `recode` stage validates and exports standardised input.

Outputs are written to `results/`, with clinical outputs under `results/clinical_cox/`. Set `CEBP_DATA_DIR` and `CEBP_RESULTS_DIR` to use other directories. `CEBP_CLINICAL_CONFIG` accepts an alternative clinical configuration file.

## Methods and outputs

- Cox models: HRs, 95% confidence intervals and two-sided P values.
- Screening: Benjamini-Hochberg adjustment across exposure contrasts.
- Competing risks: cause-specific Cox models, Fine-Gray sHRs and Aalen-Johansen estimates.
- Absolute risk, bootstrap uncertainty, concordance and proportional-hazards diagnostics.
- Multiple imputation, Rubin pooling, selection weighting and censoring weighting.
- Latent class comparison, posterior assignment and bootstrap stability.

The clinical runner reports covariate omissions and model warnings. Imputation and weighting are configured separately through the general analysis functions. Join latent class assignments to survival data by unique ID.

## Tests

```sh
Rscript tests/test_clinical.R
```

Integration tests use the public `survival::lung` dataset to check model estimates, workflow outputs and input validation.

## Data access

Participant-level study data are not included. Input data require the appropriate access authorisation.

## License

A reuse license has not yet been assigned. Contact the repository owner regarding reuse permissions.
