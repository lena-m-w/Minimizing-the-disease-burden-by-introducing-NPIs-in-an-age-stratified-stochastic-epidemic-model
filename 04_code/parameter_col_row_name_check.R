# check if parameters have the right col and row name 
# check U (interventions) and W scaling matrix


# ── Check u ───────────────────────────────────────────────────────
if (is.null(names(u))) {
  stop("u must be a named vector, e.g. c('0-19' = 0.9, '20-64' = 0.3, '65+' = 0.9)")
}

if (any(u < 0) || any(u > 1)) {
  stop("All entries of u must be in [0, 1]")
}

n <- length(u)
age_names <- names(u)

# ── Without W ─────────────────────────────────────────────────────
if (is.null(W)) {
  s <- sqrt(1 - u)
  U <- outer(s, s)
  rownames(U) <- age_names
  colnames(U) <- age_names
  return(U)
}

# ── Check W ───────────────────────────────────────────────────────
if (is.null(rownames(W)) || is.null(colnames(W))) {
  stop("W must have row and column names matching the names of u")
}

if (!all(rownames(W) == colnames(W))) {
  stop("W row and column names must be in the same order")
}

if (!all(age_names == rownames(W))) {
  stop(paste0(
    "Names of u and W do not match or are in different order.\n",
    "  names(u): ", paste(age_names, collapse = ", "), "\n",
    "  rownames(W): ", paste(rownames(W), collapse = ", ")
  ))
}

if (max(abs(W - t(W))) > tol) {
  stop("W must be symmetric (W_ab = W_ba)")
}

if (any(W < 0) || any(W > 1)) {
  stop("All entries of W must be in [0, 1]")
}