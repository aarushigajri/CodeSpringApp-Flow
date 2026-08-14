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
assert_true(sarek_manifest_validation_kind("No manifest has been validated yet.") == "info", "Pending validation should use the information style.")
assert_true(sarek_manifest_validation_kind("VALID\n- Ready") == "success", "Valid output should use the success style.")
assert_true(sarek_manifest_validation_kind("VALID\n- Ready\n\nWARNINGS\n- Review") == "warning", "Warnings should override the success style.")
assert_true(sarek_manifest_validation_kind("ERRORS\n- Fix this") == "error", "Validation errors should use the error style.")

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

sample_review <- sarek_sample_review_table(table)
assert_true(NROW(sample_review) == 1L, "File rows were not condensed into one sample-level review row.")
assert_true(sample_review$FASTQ_pairing[[1]] == "Complete: 1 lane", "Complete paired FASTQs were not recognized.")
assert_true(startsWith(sarek_fastq_pairing_status(table[1, , drop = FALSE]), "Needs attention:"), "A missing FASTQ mate was not surfaced before validation.")
assert_true(sarek_recommend_analysis_mode(table) == "tumor_only", "Tumor FASTQs did not recommend tumor-only mode.")

sample_key <- sarek_sample_key(table$patient_id[[1]], table$sample_id[[1]])
sample_updated <- sarek_apply_sample_update(
  table,
  sample_key = sample_key,
  include = TRUE,
  patient_id = "patient",
  sample_id = "patient_tumor",
  role = "tumor",
  matched_normal_id = "",
  processing_state = "unmapped"
)
assert_true(all(sample_updated$sample_id == "patient_tumor"), "A sample-level edit was not applied to every associated file.")

file_updated <- sarek_apply_file_pairing_update(
  table,
  path = table$path[[1]],
  lane = "L999",
  read = 2L
)
assert_true(file_updated$lane[[1]] == "L999" && file_updated$read[[1]] == 2L, "A FASTQ lane/read correction was not applied.")

normal_r1 <- file.path(test_root, "patient_N_L001_R1.fastq.gz")
normal_r2 <- file.path(test_root, "patient_N_L001_R2.fastq.gz")
file.create(normal_r1, normal_r2)
matched_table <- sarek_build_discovery_table(
  c(fastq_r1, fastq_r2, normal_r1, normal_r2),
  allowed_roots = test_root
)
assert_true(sarek_recommend_analysis_mode(matched_table) == "matched_tumor_normal", "Tumor-normal roles did not recommend matched mode.")

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
assert_true(grepl("sarek-review-checklist", ui_text, fixed = TRUE), "The pre-validation checklist is missing.")
assert_true(grepl("sarek_test-sample_editor", ui_text, fixed = TRUE), "The sample-level editor output is missing.")
assert_true(grepl("sarek_test-file_editor", ui_text, fixed = TRUE), "The FASTQ lane/read correction output is missing.")
assert_true(grepl("Field definitions and rules", ui_text, fixed = TRUE), "The manifest field guide is missing.")
assert_true(grepl("Results root", ui_text, fixed = TRUE) && grepl("Work root", ui_text, fixed = TRUE), "Storage fields are not explained.")
assert_true(grepl("click any sample row once", ui_text, ignore.case = TRUE), "The table does not clearly explain how to edit a sample.")
assert_true(grepl("sarek_test-validation", ui_text, fixed = TRUE), "The visible validation panel output is missing.")

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
    editor_html <- paste(as.character(output$sample_editor), collapse = "\n")
    assert_true(grepl("edit_include", editor_html, fixed = TRUE), "The selected-sample editor does not expose an include checkbox.")

    session$setInputs(sample_review_table_rows_selected = 1L)
    session$flushReact()
    expected_key <- sarek_sample_key(confirmation_state()$patient_id[[1]], confirmation_state()$sample_id[[1]])
    assert_true(selected_sample_state() == expected_key, "A single table-row selection did not select the sample for editing.")

    session$setInputs(confirm = 1L)
    session$flushReact()
    manifest <- confirmed_manifest()
    assert_true(!is.null(manifest), "Valid reviewed inputs did not produce a confirmed manifest.")
    validation_html <- paste(as.character(output$validation), collapse = "\n")
    assert_true(grepl("sarek-validation-panel-success", validation_html, fixed = TRUE), "Successful validation did not produce a green success panel.")
    download_html <- paste(as.character(output$download_ui), collapse = "\n")
    assert_true(grepl("validation_reviewed", download_html, fixed = TRUE), "Validation acknowledgement is missing before download.")
    assert_true(!grepl("download_manifest", download_html, fixed = TRUE), "Download was enabled before validation acknowledgement.")
    session$setInputs(validation_reviewed = TRUE)
    session$flushReact()
    download_html <- paste(as.character(output$download_ui), collapse = "\n")
    assert_true(grepl("download_manifest", download_html, fixed = TRUE), "Download was not enabled after validation acknowledgement.")
    assert_true(manifest$analysis$mode == "tumor_only", "Confirmed manifest has the wrong analysis mode.")
    assert_true(manifest$manifest_id == "tumor_test", "Confirmed manifest has the wrong ID.")
    assert_true(length(manifest$patients[[1]]$samples[[1]]$files) == 2L, "Confirmed manifest lost a FASTQ mate.")

    session$setInputs(preset = "changed_after_confirmation")
    session$flushReact()
    assert_true(is.null(confirmed_manifest()), "Changing settings should invalidate the prior confirmation.")
  }
)

cat("Sarek Shiny manifest tests: PASS\n")
