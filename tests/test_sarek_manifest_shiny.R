args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_manifest_shiny.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

if (!requireNamespace("shiny", quietly = TRUE)) stop("The shiny package is required for this test.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("The jsonlite package is required for this test.")

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_manifest_shiny.R"), local = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

paths <- sarek_parse_path_input(" /tmp/a.fastq.gz\n\n/tmp/b.fastq.gz\n/tmp/a.fastq.gz ")
assert_true(identical(paths, c("/tmp/a.fastq.gz", "/tmp/b.fastq.gz")), "Path input was not trimmed and deduplicated.")
assert_true(isTRUE(sarek_parse_include_value("yes")), "Truthy include value was not recognized.")
assert_true(identical(sarek_parse_include_value("0"), FALSE), "False include value was not recognized.")
assert_true(is.na(sarek_parse_include_value("maybe")), "Ambiguous include value should remain invalid.")

widths <- sarek_manifest_column_widths()
assert_true(all(c("patient_id", "sample_id", "path", "warning") %in% names(widths)), "Manifest column widths are incomplete.")
assert_true(identical(widths[["path"]], "360px"), "The input path column should remain readable.")
assert_true(sarek_manifest_status_kind("Enter input paths.") == "info", "Instruction status should use the information style.")
assert_true(sarek_manifest_status_kind("ACTION REQUIRED: Review fields.") == "review", "Review status should use the attention style.")
assert_true(sarek_manifest_status_kind("ERROR: Missing path.") == "error", "Error status should use the error style.")

test_root <- tempfile("sarek-shiny-test-")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

fastq_r1 <- file.path(test_root, "patient_T_L001_R1.fastq.gz")
fastq_r2 <- file.path(test_root, "patient_T_L001_R2.fastq.gz")
file.create(fastq_r1, fastq_r2)

table <- sarek_build_discovery_table(test_root, allowed_roots = test_root)
assert_true(NROW(table) == 2L, "Synthetic paired FASTQs were not discovered.")
summary <- sarek_confirmation_summary(table)
assert_true(summary$Value[summary$Measure == "Files"] == "2", "Discovery summary file count is incorrect.")
assert_true(summary$Value[summary$Measure == "Samples"] == "1", "Discovery summary sample count is incorrect.")

patient_col <- match("patient_id", names(table)) - 1L
edited <- sarek_apply_confirmation_edit(table, list(row = 1L, col = patient_col, value = "patient_renamed"))
assert_true(edited$patient_id[[1]] == "patient_renamed", "Editable confirmation field was not updated.")

path_col <- match("path", names(table)) - 1L
readonly <- sarek_apply_confirmation_edit(table, list(row = 1L, col = path_col, value = "/tmp/replaced.fastq.gz"))
assert_true(identical(readonly$path, table$path), "Read-only file path was unexpectedly editable.")

include_col <- match("include", names(table)) - 1L
invalid_include <- sarek_apply_confirmation_edit(table, list(row = 1L, col = include_col, value = "maybe"))
assert_true(is.na(invalid_include$include[[1]]), "Ambiguous include edit should remain invalid for validation.")
include_validation <- sarek_validate_confirmation_table(invalid_include)
assert_true(!include_validation$valid, "Invalid include edit unexpectedly passed validation.")

ui_text <- paste(as.character(sarek_manifest_ui("sarek_test")), collapse = "\n")
assert_true(grepl("sarek_test-paths", ui_text, fixed = TRUE), "Shiny module inputs were not namespaced.")
assert_true(grepl("corroborate", ui_text, ignore.case = TRUE), "UI does not explain manual corroboration.")
assert_true(grepl("sarek-confirmation-table", ui_text, fixed = TRUE), "Manifest table is missing its horizontal-scroll container.")
assert_true(grepl("sarek-status-banner", ui_text, fixed = TRUE), "Manifest status is missing its visible banner style.")

results_root <- file.path(test_root, "results")
work_root <- file.path(test_root, "work")

shiny::testServer(
  sarek_manifest_server,
  args = list(
    default_results_root = results_root,
    default_work_root = work_root,
    created_by = "test_user",
    allowed_input_roots = test_root,
    allowed_results_roots = test_root,
    allowed_work_roots = test_root,
    max_files = 10L
  ),
  {
    session$setInputs(
      paths = paste(fastq_r1, fastq_r2, sep = "\n"),
      recursive = FALSE,
      manifest_id = "tumor_test",
      assay_type = "WGS",
      analysis_mode = "tumor_only",
      preset = "core",
      results_root = results_root,
      work_root = work_root,
      species = "human",
      assembly = "GRCh38",
      sarek_genome = "GATK.GRCh38"
    )
    session$setInputs(discover = 1L)
    session$flushReact()
    assert_true(NROW(confirmation_state()) == 2L, "Module discovery did not populate the confirmation table.")
    assert_true(is.null(confirmed_manifest()), "Discovery should not automatically confirm a manifest.")

    session$setInputs(confirm = 1L)
    session$flushReact()
    manifest <- confirmed_manifest()
    assert_true(!is.null(manifest), "Valid reviewed inputs did not produce a confirmed manifest.")
    assert_true(manifest$analysis$mode == "tumor_only", "Confirmed manifest has the wrong analysis mode.")
    assert_true(manifest$manifest_id == "tumor_test", "Confirmed manifest has the wrong ID.")
    assert_true(length(manifest$patients[[1]]$samples[[1]]$files) == 2L, "Confirmed manifest lost a FASTQ mate.")

    session$setInputs(preset = "changed_after_confirmation")
    session$flushReact()
    assert_true(is.null(confirmed_manifest()), "Changing settings should invalidate the prior confirmation.")
  }
)

cat("Sarek Shiny manifest tests: PASS\n")
