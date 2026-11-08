## Validate on the same inputs as the bugphyzzExports export: drop the 2020
## Bergey's scrape (bugphyzzExports ADR 0002, superseded by Bergeys_2026) and
## the rows in sources/sheet_overrides.tsv (ADR 0005). Rules are read from a
## bugphyzzExports clone, BUGPHYZZEXPORTS_DIR (default ~/git/bugphyzzExports).
exports_dir <- path.expand(Sys.getenv("BUGPHYZZEXPORTS_DIR", "~/git/bugphyzzExports"))
source(file.path(exports_dir, "inst", "scripts", "sheet_overrides.R"))
export_inputs <- function(df, phys_name, log = message) {
    n <- sum(df$Attribute_source == "Bergey's Manual", na.rm = TRUE)
    df <- df[is.na(df$Attribute_source) | df$Attribute_source != "Bergey's Manual", ]
    log(paste0(phys_name, ": dropped ", n, " Bergey's Manual (2020) rows"))
    ov <- read_sheet_overrides(file.path(exports_dir, "sources", "sheet_overrides.tsv"))
    ov <- ov[ov$physiology == phys_name, ]
    l <- stats::setNames(list(df), phys_name)
    apply_sheet_overrides(l, ov, log = log)[[1]]
}
