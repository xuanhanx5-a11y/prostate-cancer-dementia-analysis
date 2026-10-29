# General calculation interface
Configure variable definitions, eligibility rules, time origin and covariate specifications before analysis.

Input analysis.csv uses unique id, positive time and status (0=censored, 1=dementia, 2=death before dementia). Time is elapsed from the explicitly defined analysis origin. Unknown dates/outcomes must not be converted to censoring. Covariate factors and formula terms are configurable; complete-case omission is not automatic.

methods.R supplies sequential_flow (first applicable exclusion with overlapping source flags), baseline_compare (available denominators and level-specific SMD), score_from_indicators (fixed scaling reusable across datasets), screen_cox (contrast-level BH FDR, not indicator-level omnibus FDR), cause-specific dementia/death Cox, Fine–Gray, Aalen–Johansen, net_risk, paired-cause cif_risk, bootstrap_stat, incremental model C/AIC, nested LR, Schoenfeld tests, period Cox and landmark sensitivity.

Absolute risks require an explicit newdata reference or averaging population and horizon in the input time unit. net_risk is 1-S from a cause-specific Cox model; it is not a competing-risk CIF. cif_risk integrates both cause-specific hazards and supports unstratified models. Bootstrap confidence intervals require resampling the entire desired fitting procedure; failed replicates are reported. Concordance is apparent Harrell concordance.

multiple_impute accepts imputed columns and predictor columns, never imputes id/time/status, adds event/death/log-time/cumulative-hazard auxiliaries, and preserves observed values. Inspect mice diagnostics. pool_cox uses Rubin pooling. Recalculate derived scores and class assignments under the declared frozen/refitted strategy in every completed dataset before pooling. Do not average individual HRs. Eligibility thresholds belong in caller rules, not in these functions.

selection_ipw fits an inclusion model in the eligible source registry; configurable trimming and weight diagnostics are returned. Apply weights to selected records by ID and use robust Cox variance. ipcw_cox requires a separate observed loss_to_followup flag, estimates conditional uncensoring probabilities and uses time-updated weights at event times with clustered robust variance. Administrative censoring is not labelled loss. Positivity truncation is explicit. Informative-missingness sensitivity uses explicitly specified departures from the imputation assumptions.

LCA uses a separate module. Output filenames identify the analysis and estimand.


# LCA input contract
Categorical latent class analysis with class-comparison metrics, posterior summaries and bootstrap stability.

Run `Rscript lca/run_lca.R clinical` or `Rscript lca/run_lca.R ukb` from code. Supply lca_clinical.csv or lca_ukb.csv in data. Optional unique `id`; remaining columns must follow this order:

Clinical: smoking, sedentary, falls, hypertension, diabetes, glucose, haemoglobin, creatinine.
UKB: alcohol, smoking, fev1, fvc, grip, hba1c, tyg, glucose, cystatin, crp.

Each indicator is 0=not adverse, 1=adverse, with no missing values. The approved indicator definitions must be applied upstream; do not infer that differently labelled categories are interchangeable. Outcomes and follow-up do not enter LCA. Missing indicator values are not imputed automatically.

Compares K=2–6 with 20 random starts. Selected K=3 clinically and K=4 in UKB follows the manuscript specification; BIC's preference is separately reported. Default bootstrap is 100 replicates with 5 starts and optimal permutation matching of profiles. Failed/iteration-limit replicates are explicitly recorded. These runtime choices are configurable analysis settings. Environment overrides: LCA_NREP, LCA_MAXITER, LCA_BOOTSTRAPS, LCA_BOOT_NREP. Low values are for smoke checks only.

Do not rename classes as low/metabolic/renal solely by numerical class index; review the conditional profiles. Individual assignment outputs remain participant-level data and should not automatically be made public.

# Clinical entry point
Run `Rscript code/main.R clinical` after setting code/clinical/cox_config.R (or CEBP_CLINICAL_CONFIG). See README for the shared id/time/status contract, factor reference levels, missingness policy and individual stages. Clinical outputs are cause-specific Cox HRs; Fine–Gray separately reports sHRs. Screening FDR includes exposure contrasts only. The validated clinical export can be supplied to the general pipeline by setting its input and covariate configuration. Join LCA assignments by unique ID, never by row order.
