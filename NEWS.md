# gwasurvivr 2.0.0

## Repair & correctness

* **Statistical validation.** Added a per-SNP oracle test suite asserting that
  gwasurvivr's fast warm-started `coxph.fit` path matches an independent
  `survival::coxph()` formula fit (coefficients, standard errors, hazard ratios
  within 1e-5; p-values within 1e-4), across covariate, no-covariate, and
  interaction models.
* Fixed analysis with no covariates (`covariates = NULL`), which previously
  errored in `rnorm()` / per-SNP fitting (#6).
* Fixed the blank "SNPs analyzed" count reported by the VCF (Michigan/Sanger)
  path (#25).
* Fixed empty output for SNP-by-covariate interaction models (#32).
* Fixed the `.snps_removed` file, which listed the wrong SNPs and reported a
  bogus MAF of 0 for variance-filtered SNPs (#36).
* Fixed "subscript out of bounds" when a grouped tibble / `dplyr` data frame was
  passed as `covariate.file`, and added clear errors when `sample.ids` do not
  match `id.column` (#19, #39).

## New features

* **Hard-call genotype support.** `michiganCoxSurv()` and `sangerCoxSurv()` gain
  a `genotype.field` argument (`"DS"`, `"GT"`, or `"HDS"`). `"GT"` converts hard
  genotype calls to ALT-allele dosage; `"HDS"` sums per-haplotype dosages.
* **Explicit effect allele.** VCF output now includes `EFFECT_ALLELE` (the ALT
  allele the hazard ratio refers to) and `OTHER_ALLELE` (REF). The `flip.dosage`
  behaviour for PLINK/GDS/IMPUTE2 is now documented precisely (default
  `TRUE` ⇒ effect allele is `A0`; `FALSE` ⇒ `A1`).
* **Left-truncated / interval survival.** All `*CoxSurv()` functions gain a
  `start.time` argument naming a per-sample entry-time column. When supplied,
  models use counting-process `Surv(start.time, time.to.event, event)` intervals
  (fit with survival's `agreg.fit`); otherwise standard right-censored survival
  is used. Validated to match `survival::coxph()` per SNP.

## Internals & robustness

* `survival::coxph.fit` is now resolved dynamically and passed to parallel
  workers, so the package no longer hard-depends on it being an exported symbol
  (robust to future `survival` releases) and never uses `:::`.
* Removed an experimental, broken "Modified" code path and assorted debugging
  cruft. Modernized `DESCRIPTION` (`Authors@R`, `BugReports`), added continuous
  integration, and refreshed documentation.
