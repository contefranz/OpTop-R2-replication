# =============================================================================
# prep_mdna.R -- build the MD&A (10-K item 7) corpus for the real-data study.
#
# Input:  Data/corpus_item7/item7_2016_2017.qs2 -- a quanteda corpus of item-7
#         texts filed in calendar 2016-2017 (14,817 docs; docvars: cik, company,
#         filing_date, period_of_report, fyear, sic, ...). Fiscal years straggle
#         (fyear 2015: 5,835; 2016: 7,474; other: 1,508), so the study is
#         parameterized by FISCAL year -- the year the MD&A actually discusses.
#
# Output: Data/MDNA/mdna_prep_<y1>_<y2>.qs2 with
#           dtm_train / dtm_ev -- the POOLED fiscal-{y1,y2} corpus, split at
#           random (stratified by fiscal year) into train and held-out sets;
#           dv_train / dv_ev (docvars incl. FF12 industry),
#           report (all prep counts: dedupe, floors, OOV mass, dims).
#
# Usage:
#   Rscript Code/prep_mdna.R [y1=2015] [y2=2016] [w_mindoc=5] [len_floor=200]
#                            [train_frac=0.7] [seed=1970] [sample_n=0] [workers=4]
#   sample_n > 0 subsamples each fiscal year (smoke runs, e.g. sample_n=200);
#   workers drives NLPstudio::tokenize_corpus() (parallel tokenization).
# =============================================================================

suppressMessages({
  library(quanteda); library(data.table); library(Matrix); library(qs2)
})
source(here::here("Code", "R", "source_all.R"))

# --- tiny key=value CLI (this driver has no profile semantics) ---------------
.args <- commandArgs(trailingOnly = TRUE)
P <- list(y1 = 2015L, y2 = 2016L, w_mindoc = 5L, len_floor = 200L,
          train_frac = 0.7, seed = 1970L, sample_n = 0L, workers = 4L)
for (a in .args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
  if (length(kv) != 2L || !kv[1L] %in% names(P))
    stop("unknown argument '", a, "'. Valid: ",
         paste(names(P), collapse = ", "), call. = FALSE)
  P[[kv[1L]]] <- if (kv[1L] == "train_frac") as.numeric(kv[2L]) else as.integer(kv[2L])
}
stopifnot(P$y1 < P$y2, P$train_frac > 0, P$train_frac < 1)
log_msg("=== MD&A prep: fiscal y1=%d (train/eval) vs y2=%d (temporal) ===",
        P$y1, P$y2)

# --- Fama-French 12 industries from 4-digit SIC (standard definitions) -------
ff12_from_sic <- function(sic) {
  s <- suppressWarnings(as.integer(sic))
  inr <- function(a, b) !is.na(s) & s >= a & s <= b
  out <- rep("Other", length(s))
  out[inr(100, 999) | inr(2000, 2399) | inr(2700, 2749) | inr(2770, 2799) |
      inr(3100, 3199) | inr(3940, 3989)] <- "NoDur"
  out[inr(2500, 2519) | inr(2590, 2599) | inr(3630, 3659) | inr(3710, 3711) |
      inr(3714, 3714) | inr(3716, 3716) | inr(3750, 3751) | inr(3792, 3792) |
      inr(3900, 3939) | inr(3990, 3999)] <- "Durbl"
  out[inr(2520, 2589) | inr(2600, 2699) | inr(2750, 2769) | inr(3000, 3099) |
      inr(3200, 3569) | inr(3580, 3629) | inr(3700, 3709) | inr(3712, 3713) |
      inr(3715, 3715) | inr(3717, 3749) | inr(3752, 3791) | inr(3793, 3799) |
      inr(3830, 3839) | inr(3860, 3899)] <- "Manuf"
  out[inr(1200, 1399) | inr(2900, 2999)] <- "Enrgy"
  out[inr(2800, 2829) | inr(2840, 2899)] <- "Chems"
  out[inr(3570, 3579) | inr(3660, 3692) | inr(3694, 3699) | inr(3810, 3829) |
      inr(7370, 7379)] <- "BusEq"
  out[inr(4800, 4899)] <- "Telcm"
  out[inr(4900, 4949)] <- "Utils"
  out[inr(5000, 5999) | inr(7200, 7299) | inr(7600, 7699)] <- "Shops"
  out[inr(2830, 2839) | inr(3693, 3693) | inr(3840, 3859) | inr(8000, 8099)] <- "Hlth"
  out[inr(6000, 6999)] <- "Money"
  out
}

# --- load, filter fiscal years, dedupe ---------------------------------------
corp_file <- proj_path("Data", "corpus_item7", "item7_2016_2017.qs2")
stopifnot(file.exists(corp_file))
corp <- qs_read(corp_file)
dv <- as.data.table(docvars(corp))
dv[, doc_id := docnames(corp)]
n_raw <- nrow(dv)

keep_fy <- dv$fyear %in% c(P$y1, P$y2)
dv <- dv[keep_fy]
corp <- corpus_subset(corp, fyear %in% c(P$y1, P$y2))
n_fy <- nrow(dv)

# one filing per (cik, period_of_report): keep the earliest filing_date
setorder(dv, cik, period_of_report, filing_date)
dup <- duplicated(dv[, .(cik, period_of_report)])
n_dup <- sum(dup)
keep_ids <- dv$doc_id[!dup]
corp <- corpus_subset(corp, docnames(corp) %in% keep_ids)
dv <- dv[!dup]

# optional smoke subsample, seeded, per fiscal year
if (P$sample_n > 0L) {
  set.seed(P$seed)
  keep_ids <- dv[, .SD[sample(.N, min(.N, P$sample_n))], by = fyear]$doc_id
  corp <- corpus_subset(corp, docnames(corp) %in% keep_ids)
  dv <- dv[doc_id %in% keep_ids]
  log_msg("smoke subsample: %d docs per fiscal year", P$sample_n)
}
dv[, ff12 := ff12_from_sic(sic)]
log_msg("filings: %d raw -> %d in fiscal {%d,%d} -> %d after dedupe (%d dropped)",
        n_raw, n_fy, P$y1, P$y2, nrow(dv), n_dup)

# --- tokenize (EDGAR_CORPUS conventions + LDA hygiene) ------------------------
# punct/symbols/numbers out, >= 4 chars, lowercase, English stopwords out.
# NLPstudio::tokenize_corpus() chunks the corpus and tokenizes in parallel
# (PSOCK, per the package's own advisory about quanteda/C++ internals under
# FORK); `...` forwards to quanteda::tokens().
log_msg("tokenizing %d documents on %d workers (the slow step)",
        nrow(dv), P$workers)
toks <- NLPstudio::tokenize_corpus(corp, ncores = P$workers,
                                   nchunks = P$workers, socket = "PSOCK",
                                   remove_punct = TRUE, remove_symbols = TRUE,
                                   remove_numbers = TRUE, remove_url = TRUE)
toks <- tokens_tolower(toks)
toks <- tokens_keep(toks, min_nchar = 4L)
toks <- tokens_remove(toks, stopwords("en"))
dfm_all <- dfm(toks)
rm(toks, corp); invisible(gc())

# --- POOLED design: both fiscal years, one random split (stratified) ----------
set.seed(P$seed)
ids_tr <- dv[, .SD[sample(.N, round(P$train_frac * .N))], by = fyear]$doc_id
ids_ev <- setdiff(dv$doc_id, ids_tr)
log_msg("pooled split: %d train / %d held-out (stratified by fiscal year)",
        length(ids_tr), length(ids_ev))

# vocabulary from the TRAINING half only (doc-freq floor); numeric strays out
dfm_tr <- dfm_subset(dfm_all, docnames(dfm_all) %in% ids_tr)
dfm_tr <- dfm_remove(dfm_tr, pattern = "^[0-9][0-9.,]*$", valuetype = "regex")
dfm_tr <- dfm_trim(dfm_tr, min_docfreq = P$w_mindoc, docfreq_type = "count")
vocab <- featnames(dfm_tr)
log_msg("training vocabulary: %d terms (doc-freq >= %d on %d training docs)",
        length(vocab), P$w_mindoc, length(ids_tr))

# project the held-out set onto the training vocabulary; record OOV mass
d_ev <- dfm_subset(dfm_all, docnames(dfm_all) %in% ids_ev)
tot_ev <- sum(d_ev)
d_ev <- dfm_match(d_ev, features = vocab)
oov_mass <- 1 - sum(d_ev) / tot_ev
log_msg("OOV mass dropped by vocabulary projection (held-out): %.4f", oov_mass)

# --- length floor AFTER projection; convert to plain sparse matrices ----------
as_dtm <- function(d) {
  m <- as(d, "CsparseMatrix")
  dimnames(m) <- list(docnames(d), featnames(d))
  m
}
floor_keep <- function(m) Matrix::rowSums(m) >= P$len_floor
dtm_train <- as_dtm(dfm_tr); k1 <- floor_keep(dtm_train)
dtm_ev    <- as_dtm(d_ev);   k2 <- floor_keep(dtm_ev)
log_msg("length floor (%d tokens): train %d->%d, held-out %d->%d",
        P$len_floor, nrow(dtm_train), sum(k1), nrow(dtm_ev), sum(k2))
dtm_train <- dtm_train[k1, , drop = FALSE]
dtm_ev    <- dtm_ev[k2, , drop = FALSE]

dv_of <- function(m) dv[match(rownames(m), doc_id)]
report <- list(
  params = P, n_raw = n_raw, n_fiscal = n_fy, n_dedupe_dropped = n_dup,
  W = length(vocab), J_train = nrow(dtm_train), J_ev = nrow(dtm_ev),
  year_mix_train = table(dv_of(dtm_train)$fyear),
  year_mix_ev = table(dv_of(dtm_ev)$fyear),
  oov_mass = oov_mass,
  len_summary_train = summary(Matrix::rowSums(dtm_train)),
  built = format(Sys.time())
)
out <- list(dtm_train = dtm_train, dtm_ev = dtm_ev,
            dv_train = dv_of(dtm_train), dv_ev = dv_of(dtm_ev),
            report = report)
f <- p_data("MDNA", sprintf("mdna_prep_%d_%d%s.qs2", P$y1, P$y2,
                            if (P$sample_n > 0L) sprintf("_n%d", P$sample_n) else ""))
qs_save(out, f)
log_msg("saved %s (train %d x %d; held-out %d)",
        basename(f), nrow(dtm_train), ncol(dtm_train), nrow(dtm_ev))
log_msg("=== MD&A prep complete ===")
