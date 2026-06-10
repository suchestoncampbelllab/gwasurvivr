# Resolvers for survival's low-level per-iteration fitters.
#
# gwasurvivr's speed depends on calling survival's internal fitters directly
# with a warm-start INIT, rather than the slow formula interface. Two fitters
# are needed, and they share the same positional signature:
#   coxph.fit() - right-censored Surv(time, event) data.
#   agreg.fit() - left-truncated / counting-process Surv(start, stop, event)
#                 data (used when start.time is supplied).
# These have historically moved between exported and internal across survival
# releases, so we resolve them dynamically (preferring the exported binding,
# falling back to the namespace) and pass the resolved closure through
# cox.params to every parallel worker. This avoids `survival:::` (which R CMD
# check / Bioconductor flag) and avoids depending on the survival namespace
# being attached on PSOCK workers.

.resolve_survival_fitter <- function(name) {
    fn <- tryCatch(getExportedValue("survival", name), error = function(e) NULL)
    if (is.null(fn)) {
        fn <- tryCatch(getFromNamespace(name, "survival"),
                       error = function(e) NULL)
    }
    if (is.null(fn) || !is.function(fn)) {
        stop("Could not locate survival::", name, "(). Please install/upgrade ",
             "the 'survival' package.")
    }
    fn
}

# Right-censored fitter (the default hot path).
.resolve_coxph_fit <- function() .resolve_survival_fitter("coxph.fit")

# Counting-process / left-truncated fitter (used when start.time is supplied).
.resolve_agreg_fit <- function() .resolve_survival_fitter("agreg.fit")
