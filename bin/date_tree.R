#!/usr/bin/env Rscript

# Time-calibrate a species tree with ape::chronos.
#
# Ported from Eco-Flow/excon (tree_subsampling branch, bin/date_tree.R), which
# adapted it from make_dated_trees.R by Tom Wenseleers (KU Leuven, 2026), which
# dated the Vespidae and Aculeata subsets of a 72-species EXCON analysis. The
# calibration ages, clade definitions and subsetting in that script were specific
# to those data; here (as in excon) the same approach is driven by a
# user-supplied calibration table so it applies to any dataset. The validation
# it performs -- treating chronos warnings as fatal and checking that the fitted
# node ages actually match the calibrations -- is kept, as is the idea of
# recording calibration provenance alongside the tree.
#
# Converts a tree whose branch lengths are substitutions per site into an
# ultrametric tree whose branch lengths are millions of years, using node ages
# supplied by the user.
#
# Arguments (positional):
#   args[1] = input tree (rooted Newick, branch lengths in substitutions/site)
#   args[2] = calibration table (TSV, see below)
#   args[3] = chronos model: "discrete", "correlated" or "relaxed"
#   args[4] = chronos lambda (rate smoothing)
#   args[5] = number of rate categories (used by the "discrete" model)
#
# The calibration table needs a header and these columns:
#   clade    label used in the report
#   tips     comma-separated tip names; their MRCA is the calibrated node
#   age_min  minimum age in Ma
#   age_max  maximum age in Ma (equal to age_min to fix the age)
#
# Naming a node by the MRCA of two tips keeps the table readable and stable:
# it does not depend on node numbering, which differs between trees.
#
# Unlike excon's version, tip names here are matched exactly (no ".clean"
# suffix stripping): that suffix comes from excon's own RENAME_FASTA step,
# which this pipeline doesn't have -- species names here are the --input
# samplesheet IDs, used as-is throughout.
#
# Outputs: SpeciesTree_dated.nwk, dating_calibrations.tsv, dating_qc.tsv

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
tree_file  <- args[1]
calib_file <- args[2]
model      <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "discrete"
lambda     <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 1
rate_cats  <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 10

cat("================================================\n")
cat("DATE TREE (ape::chronos)\n")
cat("Model:", model, " Lambda:", lambda, " Rate categories:", rate_cats, "\n")
cat("================================================\n")

tre <- read.tree(tree_file)
stopifnot("Input tree must be rooted" = is.rooted(tre))
if (is.null(tre$edge.length)) stop("Input tree has no branch lengths.")
if (any(!is.finite(tre$edge.length)) || any(tre$edge.length <= 0)) {
  stop("All input branch lengths must be finite and positive.")
}

# Support values sit in node labels and confuse chronos; they are not needed
# for dating and the original tree is published alongside anyway.
tre$node.label <- NULL

calib <- read.delim(calib_file, stringsAsFactors = FALSE)
required <- c("clade", "tips", "age_min", "age_max")
missing_cols <- setdiff(required, names(calib))
if (length(missing_cols)) {
  stop("Calibration table is missing column(s): ", paste(missing_cols, collapse = ", "))
}
if (!nrow(calib)) stop("Calibration table is empty.")

resolve_tips <- function(spec, clade) {
  wanted <- trimws(strsplit(spec, ",")[[1]])
  wanted <- wanted[nzchar(wanted)]
  if (length(wanted) < 2) {
    stop("Calibration '", clade, "' needs at least two tips to define an MRCA.")
  }
  missing <- setdiff(wanted, tre$tip.label)
  if (length(missing)) {
    stop("Calibration '", clade, "' names tips absent from the tree: ",
         paste(missing, collapse = ", "),
         "\nTree tips are: ", paste(tre$tip.label, collapse = ", "))
  }
  wanted
}

nodes <- integer(nrow(calib))
for (i in seq_len(nrow(calib))) {
  nodes[i] <- getMRCA(tre, resolve_tips(calib$tips[i], calib$clade[i]))
}
if (anyDuplicated(nodes)) {
  dup <- calib$clade[duplicated(nodes) | duplicated(nodes, fromLast = TRUE)]
  stop("These calibrations resolve to the same node: ", paste(dup, collapse = ", "))
}

cat("Calibrations:\n")
for (i in seq_len(nrow(calib))) {
  cat(sprintf("   %-28s node %-5d %.2f - %.2f Ma\n",
              calib$clade[i], nodes[i], calib$age_min[i], calib$age_max[i]))
}

calibration <- makeChronosCalib(
  tre,
  node        = nodes,
  age.min     = calib$age_min,
  age.max     = calib$age_max,
  soft.bounds = FALSE
)

fit_warnings <- character()
fit <- withCallingHandlers(
  chronos(
    tre,
    lambda      = lambda,
    model       = model,
    calibration = calibration,
    control     = chronos.control(nb.rate.cat = rate_cats,
                                  iter.max = 20000, eval.max = 20000,
                                  dual.iter.max = 50),
    quiet       = TRUE
  ),
  warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

# chronos reports non-convergence as a warning, so treat warnings as fatal
# rather than silently publishing a tree that did not converge.
if (length(fit_warnings)) {
  stop("chronos reported: ", paste(unique(fit_warnings), collapse = " | "))
}

class(fit) <- "phylo"
stopifnot("Dated tree must be ultrametric" = is.ultrametric(fit, tol = 1e-7))
stopifnot("Dated tree must be binary"      = is.binary(fit))
if (any(!is.finite(fit$edge.length)) || any(fit$edge.length <= 0)) {
  stop("Dated tree contains non-positive or non-finite branch lengths.")
}

node_age <- function(tree, node) {
  d <- node.depth.edgelength(tree)
  max(d[seq_len(Ntip(tree))]) - d[node]
}

fitted_nodes <- integer(nrow(calib))
for (i in seq_len(nrow(calib))) {
  fitted_nodes[i] <- getMRCA(fit, resolve_tips(calib$tips[i], calib$clade[i]))
}
fitted_ages <- vapply(fitted_nodes, function(n) node_age(fit, n), numeric(1))

tol <- 1e-4
within <- fitted_ages >= calib$age_min - tol & fitted_ages <= calib$age_max + tol
if (!all(within)) {
  stop("The dated tree does not reproduce these calibrations: ",
       paste(sprintf("%s (wanted %.2f-%.2f, got %.2f)",
                     calib$clade[!within], calib$age_min[!within],
                     calib$age_max[!within], fitted_ages[!within]),
             collapse = "; "))
}

write.tree(fit, "SpeciesTree_dated.nwk")

report <- data.frame(
  clade         = calib$clade,
  tips          = calib$tips,
  node          = fitted_nodes,
  age_min_Ma    = calib$age_min,
  age_max_Ma    = calib$age_max,
  fitted_age_Ma = round(fitted_ages, 6),
  within_bounds = within,
  stringsAsFactors = FALSE
)
write.table(report, "dating_calibrations.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

depths <- node.depth.edgelength(fit)
tip_depths <- depths[seq_len(Ntip(fit))]
qc <- data.frame(
  input_tree       = tree_file,
  model            = model,
  lambda           = lambda,
  rate_categories  = rate_cats,
  n_tips           = Ntip(fit),
  n_calibrations   = nrow(calib),
  root_age_Ma      = round(max(tip_depths), 6),
  min_edge_Ma      = round(min(fit$edge.length), 8),
  ultrametric      = is.ultrametric(fit, tol = 1e-7),
  binary           = is.binary(fit),
  stringsAsFactors = FALSE
)
write.table(qc, "dating_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("------------------------------------------------\n")
cat("Dated tree written: SpeciesTree_dated.nwk\n")
cat("   Tips:      ", Ntip(fit), "\n")
cat("   Root age:  ", round(max(tip_depths), 3), "Ma\n")
cat("   Min branch:", signif(min(fit$edge.length), 3), "Ma\n")
cat("================================================\n")
