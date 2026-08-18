# Shiny adapter for the pure Sarek manifest helpers in sarek_manifest.R.
# Keeping UI state here avoids adding pipeline-specific logic to app.R.

sarek_shiny_value <- function(value, default = "") {
  if (is.null(value) || !length(value) || is.na(value[[1]])) return(default)
  value
}

sarek_parse_path_input <- function(value) {
  lines <- unlist(strsplit(as.character(sarek_shiny_value(value)), "\n", fixed = TRUE), use.names = FALSE)
  paths <- trimws(lines)
  unique(paths[!is.na(paths) & nzchar(paths)])
}

sarek_parse_include_value <- function(value) {
  value <- tolower(trimws(as.character(sarek_shiny_value(value))))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  NA
}

sarek_editable_columns <- function() {
  c(
    "include",
    "patient_id",
    "sample_id",
    "role",
    "matched_normal_id",
    "input_format",
    "processing_state",
    "lane",
    "read"
  )
}

sarek_manifest_column_widths <- function() {
  c(
    include = "80px",
    patient_id = "150px",
    sample_id = "170px",
    role = "110px",
    matched_normal_id = "180px",
    input_format = "120px",
    processing_state = "170px",
    path = "360px",
    index = "320px",
    lane = "90px",
    read = "80px",
    size_bytes = "120px",
    role_confidence = "140px",
    warning = "340px"
  )
}

sarek_manifest_column_defs <- function(table) {
  widths <- sarek_manifest_column_widths()
  widths <- widths[names(widths) %in% names(table)]
  lapply(names(widths), function(column) {
    list(targets = match(column, names(table)) - 1L, width = unname(widths[[column]]))
  })
}

sarek_manifest_status_kind <- function(message) {
  message <- trimws(as.character(sarek_shiny_value(message)))
  if (startsWith(message, "ERROR:")) return("error")
  if (startsWith(message, "ACTION REQUIRED:")) return("review")
  "info"
}

sarek_sample_key <- function(patient_id, sample_id) {
  paste(as.character(patient_id), as.character(sample_id), sep = "::")
}

sarek_fastq_pairing_status <- function(rows) {
  include <- suppressWarnings(as.logical(rows$include))
  rows <- rows[which(!is.na(include) & include), , drop = FALSE]
  if (!NROW(rows)) return("Excluded")
  rows <- rows[rows$input_format == "fastq", , drop = FALSE]
  if (!NROW(rows)) return("Not applicable")

  lane_values <- trimws(as.character(rows$lane))
  lanes <- ifelse(is.na(lane_values) | !nzchar(lane_values), "unlabelled", lane_values)
  reads_by_lane <- split(rows$read, lanes)
  problems <- unlist(lapply(names(reads_by_lane), function(lane) {
    reads <- suppressWarnings(as.integer(reads_by_lane[[lane]]))
    if (length(reads) == 2L && !anyNA(reads) && identical(sort(reads), c(1L, 2L))) {
      return(character(0))
    }
    details <- character(0)
    missing <- setdiff(c(1L, 2L), reads[!is.na(reads)])
    if (length(missing)) details <- c(details, paste0("missing R", missing))
    if (anyNA(reads)) details <- c(details, "read number not detected")
    duplicated_reads <- unique(reads[!is.na(reads) & duplicated(reads)])
    if (length(duplicated_reads)) details <- c(details, paste0("duplicate R", duplicated_reads))
    if (length(reads) > 2L && !length(duplicated_reads)) details <- c(details, "more than two files")
    paste0(lane, " ", paste(details, collapse = ", "))
  }), use.names = FALSE)

  if (length(problems)) paste0("Needs attention: ", paste(problems, collapse = "; ")) else {
    paste0("Complete: ", length(reads_by_lane), " lane", if (length(reads_by_lane) == 1L) "" else "s")
  }
}

sarek_sample_review_table <- function(table) {
  if (is.null(table) || !NROW(table)) return(data.frame())
  keys <- sarek_sample_key(table$patient_id, table$sample_id)
  groups <- split(seq_len(NROW(table)), keys)
  rows <- lapply(groups, function(index) {
    sample_rows <- table[index, , drop = FALSE]
    values <- function(column, blank = "") {
      value <- unique(trimws(as.character(sample_rows[[column]])))
      value <- value[!is.na(value) & nzchar(value)]
      if (length(value)) paste(value, collapse = ", ") else blank
    }
    data.frame(
      include = all(suppressWarnings(as.logical(sample_rows$include)), na.rm = TRUE),
      patient_id = as.character(sample_rows$patient_id[[1]]),
      sample_id = as.character(sample_rows$sample_id[[1]]),
      role = values("role", "unknown"),
      matched_normal_id = values("matched_normal_id"),
      input_format = values("input_format"),
      processing_state = values("processing_state"),
      file_count = NROW(sample_rows),
      FASTQ_pairing = sarek_fastq_pairing_status(sample_rows),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sarek_recommend_analysis_mode <- function(table) {
  if (is.null(table) || !NROW(table)) return("")
  include <- suppressWarnings(as.logical(table$include))
  selected <- table[which(!is.na(include) & include), , drop = FALSE]
  if (!NROW(selected)) return("")
  formats <- unique(as.character(selected$input_format))
  roles <- unique(as.character(selected$role))
  if (all(formats %in% c("vcf", "bcf"))) return("annotation_only")
  if (length(roles) == 1L && identical(roles, "germline")) return("germline")
  if (length(roles) == 1L && identical(roles, "tumor")) return("tumor_only")
  if (all(c("tumor", "normal") %in% roles)) {
    by_patient <- split(selected$role, selected$patient_id)
    complete <- vapply(by_patient, function(x) any(x == "tumor") && any(x == "normal"), logical(1))
    if (all(complete)) return("matched_tumor_normal")
  }
  ""
}

sarek_analysis_mode_label <- function(mode) {
  labels <- c(
    germline = "Germline",
    tumor_only = "Tumor only",
    matched_tumor_normal = "Matched tumor-normal",
    annotation_only = "Annotation only"
  )
  if (mode %in% names(labels)) unname(labels[[mode]]) else "No clear recommendation"
}

sarek_apply_sample_update <- function(
  table,
  sample_key,
  include,
  patient_id,
  sample_id,
  role,
  matched_normal_id,
  processing_state
) {
  keys <- sarek_sample_key(table$patient_id, table$sample_id)
  rows <- which(keys == sample_key)
  if (!length(rows)) stop("The selected sample is no longer available. Select it again.")
  table$include[rows] <- isTRUE(include)
  table$patient_id[rows] <- trimws(as.character(patient_id))
  table$sample_id[rows] <- trimws(as.character(sample_id))
  table$role[rows] <- as.character(role)
  table$matched_normal_id[rows] <- if (identical(role, "tumor")) trimws(as.character(matched_normal_id)) else ""
  table$processing_state[rows] <- as.character(processing_state)
  table
}

sarek_apply_file_pairing_update <- function(table, path, lane, read) {
  rows <- which(as.character(table$path) == as.character(path))
  if (length(rows) != 1L) stop("The selected file is no longer available. Select it again.")
  if (!identical(as.character(table$input_format[[rows]]), "fastq")) {
    stop("Lane and read corrections apply only to FASTQ files.")
  }
  read <- suppressWarnings(as.integer(read))
  if (!is.na(read) && !read %in% c(1L, 2L)) stop("FASTQ read must be R1, R2, or not detected.")
  table$lane[[rows]] <- trimws(as.character(lane))
  table$read[[rows]] <- read
  table
}

sarek_apply_confirmation_edit <- function(table, edit) {
  if (is.null(table) || !NROW(table) || is.null(edit$row) || is.null(edit$col)) return(table)
  row <- suppressWarnings(as.integer(edit$row))
  column <- suppressWarnings(as.integer(edit$col)) + 1L
  if (
    is.na(row) || is.na(column) ||
    row < 1L || row > NROW(table) ||
    column < 1L || column > NCOL(table)
  ) {
    return(table)
  }

  column_name <- names(table)[[column]]
  if (!column_name %in% sarek_editable_columns()) return(table)
  value <- as.character(sarek_shiny_value(edit$value))

  if (identical(column_name, "include")) {
    table[[column_name]][[row]] <- sarek_parse_include_value(value)
  } else if (identical(column_name, "read")) {
    value <- trimws(value)
    table[[column_name]][[row]] <- if (nzchar(value)) suppressWarnings(as.integer(value)) else NA_integer_
  } else {
    table[[column_name]][[row]] <- trimws(value)
  }
  table
}

sarek_confirmation_summary <- function(table) {
  if (is.null(table) || !NROW(table)) {
    return(data.frame(
      Measure = c("Files", "Included files", "Patients", "Samples"),
      Value = c(0L, 0L, 0L, 0L),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  include <- suppressWarnings(as.logical(table$include))
  selected <- table[which(!is.na(include) & include), , drop = FALSE]
  sample_keys <- if (NROW(selected)) {
    unique(paste(selected$patient_id, selected$sample_id, sep = "::"))
  } else {
    character(0)
  }
  formats <- if (NROW(selected)) paste(sort(unique(selected$input_format)), collapse = ", ") else ""
  roles <- if (NROW(selected)) paste(sort(unique(selected$role)), collapse = ", ") else ""
  data.frame(
    Measure = c("Files", "Included files", "Patients", "Samples", "Formats", "Roles"),
    Value = c(
      as.character(NROW(table)),
      as.character(NROW(selected)),
      as.character(length(unique(selected$patient_id))),
      as.character(length(sample_keys)),
      formats,
      roles
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sarek_manifest_validation_text <- function(validation, mode_errors = character(0)) {
  errors <- unique(c(sarek_shiny_value(validation$errors, character(0)), mode_errors))
  warnings <- unique(sarek_shiny_value(validation$warnings, character(0)))
  sections <- character(0)
  if (length(errors)) {
    sections <- c(sections, "ERRORS", paste0("- ", errors))
  } else {
    sections <- c(sections, "VALID", "- The current confirmation table and analysis mode are valid.")
  }
  if (length(warnings)) {
    sections <- c(sections, "", "WARNINGS", paste0("- ", warnings))
  }
  paste(sections, collapse = "\n")
}

sarek_manifest_validation_kind <- function(message) {
  message <- as.character(sarek_shiny_value(message))
  if (grepl("(^|\\n)ERRORS($|\\n)", message, perl = TRUE)) return("error")
  if (grepl("(^|\\n)WARNINGS($|\\n)", message, perl = TRUE)) return("warning")
  if (grepl("(^|\\n)(VALID|CONFIRMED)($|\\n)", message, perl = TRUE)) return("success")
  "info"
}

sarek_manifest_validation_heading <- function(kind) {
  switch(
    kind,
    error = "Validation failed — corrections required",
    warning = "Validation passed with warnings — review required",
    success = "Validation passed",
    "Validation has not run"
  )
}

sarek_manifest_field_guide <- function() {
  shiny::tags$details(
    class = "sarek-field-guide",
    open = "open",
    shiny::tags$summary("Field definitions and rules"),
    shiny::fluidRow(
      shiny::column(
        4,
        shiny::h5("Sample identity"),
        shiny::tags$dl(
          shiny::tags$dt("Include"),
          shiny::tags$dd("Checkbox controlling whether the entire sample and all its files enter the manifest."),
          shiny::tags$dt("Patient ID"),
          shiny::tags$dd("Groups related samples. A matched tumor and normal must use the same patient ID."),
          shiny::tags$dt("Sample ID"),
          shiny::tags$dd("Identifies one biological sample within a patient. It must be unique within that patient."),
          shiny::tags$dt("Role"),
          shiny::tags$dd("Germline, tumor, normal, or unknown. Filename-based role guesses must be reviewed."),
          shiny::tags$dt("Matched normal"),
          shiny::tags$dd("For a tumor sample, select an included normal sample belonging to the same patient.")
        )
      ),
      shiny::column(
        4,
        shiny::h5("Input files"),
        shiny::tags$dl(
          shiny::tags$dt("Input format"),
          shiny::tags$dd("Detected from the file extension and read-only: FASTQ, uBAM, BAM, CRAM, VCF, or BCF."),
          shiny::tags$dt("Processing state"),
          shiny::tags$dd("How far the input has already been processed, from unmapped through analysis-ready or variant calls."),
          shiny::tags$dt("FASTQ lane and read"),
          shiny::tags$dd("Every included sample lane must contain exactly one R1 and one R2. Correct detection only when the filename was interpreted incorrectly."),
          shiny::tags$dt("Index"),
          shiny::tags$dd("Detected companion index for BAM, CRAM, VCF, or BCF. Missing indexes are reported before pipeline submission.")
        )
      ),
      shiny::column(
        4,
        shiny::h5("Analysis and storage"),
        shiny::tags$dl(
          shiny::tags$dt("Analysis mode"),
          shiny::tags$dd("Must agree with the included roles: all germline, all tumor, or matched tumor-normal."),
          shiny::tags$dt("Results root"),
          shiny::tags$dd("Permanent destination for final pipeline outputs and reports. It must be an absolute writable server path."),
          shiny::tags$dt("Work root"),
          shiny::tags$dd("High-capacity temporary space used by Nextflow for intermediate files. It can be much larger than the final results and should not use a space-limited home directory."),
          shiny::tags$dt("Manifest ID"),
          shiny::tags$dd("Short name used to identify this requested analysis and its exported JSON manifest.")
        )
      )
    )
  )
}

sarek_manifest_table_output <- function(output_id) {
  if (requireNamespace("DT", quietly = TRUE)) {
    DT::DTOutput(output_id)
  } else {
    shiny::tableOutput(output_id)
  }
}

sarek_manifest_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$style(shiny::HTML(
      ".sarek-confirmation-table { width:100%; max-width:100%; overflow-x:auto; padding-bottom:8px; }
       .sarek-confirmation-table .dataTables_wrapper { overflow-x:visible; }
       .sarek-confirmation-table .dataTables_scroll { width:100%; overflow-x:auto; }
       .sarek-confirmation-table .dataTables_scrollHeadInner,
       .sarek-confirmation-table table.dataTable {
         width:max-content !important;
         min-width:2550px !important;
         table-layout:auto !important;
       }
       .sarek-confirmation-table table.dataTable thead th,
       .sarek-confirmation-table table.dataTable tbody td {
         min-width:110px !important;
         max-width:none !important;
         white-space:nowrap !important;
       }
       .sarek-confirmation-table table.dataTable thead th:nth-child(8),
       .sarek-confirmation-table table.dataTable tbody td:nth-child(8),
       .sarek-confirmation-table table.dataTable thead th:nth-child(9),
       .sarek-confirmation-table table.dataTable tbody td:nth-child(9),
       .sarek-confirmation-table table.dataTable thead th:nth-child(14),
       .sarek-confirmation-table table.dataTable tbody td:nth-child(14) {
         min-width:320px !important;
       }
       .sarek-status-banner {
         display:flex;
         align-items:flex-start;
         gap:10px;
         width:100%;
         margin:8px 0 16px 0;
         padding:13px 15px;
         border:1px solid #b8d4f2;
         border-left:6px solid #1769aa;
         border-radius:8px;
         background:#eef6ff;
         color:#17324d;
         font-size:14px;
         font-weight:650;
         line-height:1.45;
         white-space:normal;
         overflow-wrap:anywhere;
         box-shadow:0 2px 7px rgba(23,50,77,.08);
       }
       .sarek-status-banner::before {
         content:'i';
         flex:0 0 22px;
         width:22px;
         height:22px;
         border-radius:50%;
         background:#1769aa;
         color:#fff;
         text-align:center;
         font-weight:800;
         line-height:22px;
       }
       .sarek-status-banner-review {
         border-color:#e7c66a;
         border-left-color:#b7791f;
         background:#fff8df;
         color:#5c3b08;
       }
       .sarek-status-banner-review::before {
         content:'!';
         background:#b7791f;
       }
       .sarek-status-banner-error {
         border-color:#efb0b0;
         border-left-color:#b42318;
         background:#fff1f0;
         color:#7a1c16;
       }
       .sarek-status-banner-error::before {
         content:'!';
         background:#b42318;
       }
       .sarek-review-checklist {
         margin:10px 0 16px 0;
         padding:12px 14px;
         border:1px solid #d7e0ea;
         border-radius:8px;
         background:#f8fafc;
       }
       .sarek-review-item { margin:5px 0; line-height:1.4; }
       .sarek-review-item-pass { color:#276749; }
       .sarek-review-item-review { color:#8a5700; font-weight:650; }
       .sarek-sample-editor {
         margin:8px 0 16px 0;
         padding:14px;
         border:1px solid #d7e0ea;
         border-radius:8px;
         background:#fff;
       }
       .sarek-sample-editor h4 { margin-top:0; }
       .sarek-sample-review-table { width:100%; max-width:100%; overflow-x:auto; }
       .sarek-sample-review-table table.dataTable {
         width:max-content !important;
         min-width:1450px !important;
         table-layout:auto !important;
       }
       .sarek-sample-review-table table.dataTable th,
       .sarek-sample-review-table table.dataTable td {
         min-width:120px !important;
         max-width:none !important;
         white-space:nowrap !important;
       }
       .sarek-edit-instruction {
         margin:10px 0;
         padding:10px 12px;
         border-left:5px solid #1769aa;
         border-radius:6px;
         background:#eef6ff;
         color:#17324d;
         font-weight:650;
       }
       .sarek-field-guide {
         margin:12px 0 18px 0;
         padding:12px 14px;
         border:1px solid #cfd9e5;
         border-radius:8px;
         background:#fbfcfe;
       }
       .sarek-field-guide summary { cursor:pointer; font-weight:700; }
       .sarek-field-guide h5 { margin-top:16px; font-weight:700; }
       .sarek-field-guide dt { margin-top:8px; color:#17324d; }
       .sarek-field-guide dd { margin-left:0; color:#536273; line-height:1.4; }
       .sarek-validation-panel {
         display:flex;
         align-items:flex-start;
         gap:12px;
         margin:10px 0 16px 0;
         padding:15px 17px;
         border:1px solid #b8d4f2;
         border-left:7px solid #1769aa;
         border-radius:8px;
         background:#eef6ff;
         color:#17324d;
         box-shadow:0 2px 8px rgba(23,50,77,.08);
       }
       .sarek-validation-panel::before {
         content:'i';
         flex:0 0 24px;
         width:24px;
         height:24px;
         border-radius:50%;
         background:#1769aa;
         color:#fff;
         text-align:center;
         font-weight:800;
         line-height:24px;
       }
       .sarek-validation-panel h4 { margin:1px 0 7px 0; font-weight:750; }
       .sarek-validation-panel pre {
         margin:0;
         padding:0;
         border:0;
         background:transparent;
         color:inherit;
         font-family:inherit;
         font-size:14px;
         line-height:1.5;
         white-space:pre-wrap;
         overflow-wrap:anywhere;
       }
       .sarek-validation-panel-success {
         border-color:#9fd4b2;
         border-left-color:#2f855a;
         background:#edf9f1;
         color:#22543d;
       }
       .sarek-validation-panel-success::before { content:'✓'; background:#2f855a; }
       .sarek-validation-panel-warning {
         border-color:#e7c66a;
         border-left-color:#b7791f;
         background:#fff8df;
         color:#5c3b08;
       }
       .sarek-validation-panel-warning::before { content:'!'; background:#b7791f; }
       .sarek-validation-panel-error {
         border-color:#efb0b0;
         border-left-color:#b42318;
         background:#fff1f0;
         color:#7a1c16;
       }
       .sarek-validation-panel-error::before { content:'!'; background:#b42318; }
       .sarek-storage-note {
         margin:-8px 0 12px 0;
         color:#536273;
         font-size:13px;
         line-height:1.4;
       }"
    )),
    shiny::div(
      class = "progress-header-row",
      shiny::div(
        shiny::h3("Build a Sarek input manifest"),
        shiny::tags$p(
          class = "muted",
          "Discover sequencing files, review every inferred field, and export a confirmed JSON manifest."
        )
      )
    ),
    sarek_manifest_field_guide(),
    shiny::tags$div(
      class = "read-source-note",
      shiny::tags$strong("Step 1: discover files"),
      shiny::tags$p(
        "Enter one readable absolute server path per line. A path may point to a file or folder."
      ),
      shiny::tags$p(
        class = "muted small-note",
        "Supported inputs: paired FASTQ, uBAM, BAM, CRAM, VCF, and BCF."
      )
    ),
    shiny::fluidRow(
      shiny::column(
        8,
        shiny::textAreaInput(
          ns("paths"),
          "Input files or folders",
          value = "",
          rows = 8,
          placeholder = "/absolute/path/sample_T_L001_R1.fastq.gz\n/absolute/path/sample_T_L001_R2.fastq.gz"
        ),
        shiny::checkboxInput(ns("recursive"), "Search inside subfolders", value = FALSE),
        shiny::actionButton(ns("discover"), "Discover inputs", class = "btn-primary")
      ),
      shiny::column(
        4,
        shiny::uiOutput(ns("status"))
      )
    ),
    shiny::tags$hr(),
    shiny::tags$div(
      class = "read-source-note",
      shiny::tags$strong("Step 2: review and correct samples"),
      shiny::tags$p(
        "The searchable table uses one row per biological sample. Select a row once to edit it below; changes apply consistently to every associated file."
      ),
      shiny::tags$p(
        class = "muted small-note",
        "File format, paths, indexes, sizes, confidence, and warnings remain read-only. Lane/read detection can be corrected inside file details."
      )
    ),
    shiny::tableOutput(ns("summary")),
    shiny::uiOutput(ns("readiness")),
    shiny::div(
      class = "sarek-edit-instruction",
      "Editable: click any sample row once, then use the checkbox, text fields, and dropdowns in the selected-sample editor below."
    ),
    shiny::div(
      class = "sarek-sample-review-table",
      sarek_manifest_table_output(ns("sample_review_table"))
    ),
    shiny::uiOutput(ns("sample_editor")),
    shiny::tags$details(
      shiny::tags$summary("Show files for the selected sample"),
      shiny::br(),
      shiny::uiOutput(ns("file_editor")),
      shiny::div(
        class = "sarek-confirmation-table",
        sarek_manifest_table_output(ns("file_detail_table"))
      )
    ),
    shiny::tags$hr(),
    shiny::tags$div(
      class = "read-source-note",
      shiny::tags$strong("Step 3: describe the intended analysis"),
      shiny::tags$p(
        "These settings describe the confirmed request. They do not submit nf-core/sarek."
      )
    ),
    shiny::fluidRow(
      shiny::column(
        4,
        shiny::textInput(ns("manifest_id"), "Manifest ID", value = "sarek_analysis"),
        shiny::selectInput(
          ns("assay_type"),
          "Assay",
          choices = c("Whole-genome sequencing" = "WGS", "Whole-exome sequencing" = "WES",
                      "Targeted panel" = "targeted", "Annotation only" = "annotation_only"),
          selected = "WGS"
        ),
        shiny::selectInput(
          ns("analysis_mode"),
          "Analysis mode",
          choices = c("Germline" = "germline", "Tumor only" = "tumor_only",
                      "Matched tumor-normal" = "matched_tumor_normal",
                      "Annotation only" = "annotation_only"),
          selected = "germline"
        ),
        shiny::textInput(ns("preset"), "Analysis preset", value = "core")
      ),
      shiny::column(
        4,
        shiny::textInput(
          ns("results_root"),
          "Results root — permanent outputs",
          value = "",
          placeholder = "/absolute/path/results/sarek/user"
        ),
        shiny::tags$p(
          class = "sarek-storage-note",
          "Final results, reports, and deliverables will be stored here."
        ),
        shiny::textInput(
          ns("work_root"),
          "Work root — temporary Nextflow files",
          value = "",
          placeholder = "/high-capacity/path/work/sarek/user"
        ),
        shiny::tags$p(
          class = "sarek-storage-note",
          "Intermediate files can be much larger than the final results. Use high-capacity storage, not a space-limited home directory."
        ),
        shiny::tags$p(
          class = "muted small-note",
          "Both values must be absolute server paths. Manifest creation records them but does not create folders or start Sarek."
        )
      ),
      shiny::column(
        4,
        shiny::tags$details(
          shiny::tags$summary("Reference settings"),
          shiny::br(),
          shiny::textInput(ns("species"), "Species", value = "human"),
          shiny::textInput(ns("assembly"), "Assembly", value = "GRCh38"),
          shiny::textInput(ns("sarek_genome"), "nf-core/sarek genome key", value = "GATK.GRCh38")
        ),
        shiny::br(),
        shiny::actionButton(ns("confirm"), "Validate and confirm manifest", class = "btn-primary"),
        shiny::br(),
        shiny::br(),
        shiny::uiOutput(ns("download_ui"))
      )
    ),
    shiny::h4("Validation"),
    shiny::uiOutput(ns("validation")),
    shiny::tags$details(
      shiny::tags$summary("Confirmed JSON preview"),
      shiny::br(),
      shiny::verbatimTextOutput(ns("manifest_preview"))
    )
  )
}

sarek_manifest_server <- function(
  id,
  default_results_root,
  default_work_root,
  created_by = "unknown",
  allowed_input_roots = NULL,
  allowed_results_roots = NULL,
  allowed_work_roots = NULL,
  max_files = 5000L
) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    confirmation_state <- shiny::reactiveVal(NULL)
    confirmed_manifest <- shiny::reactiveVal(NULL)
    status_state <- shiny::reactiveVal(
      "Enter one or more server paths, then select Discover inputs."
    )
    validation_state <- shiny::reactiveVal(
      "No manifest has been validated yet."
    )
    selected_sample_state <- shiny::reactiveVal("")
    validation_version <- shiny::reactiveVal(0L)
    reviewed_validation_version <- shiny::reactiveVal(-1L)

    shiny::observe({
      shiny::updateTextInput(session, "results_root", value = default_results_root)
      shiny::updateTextInput(session, "work_root", value = default_work_root)
    })

    invalidate_confirmation <- function() {
      confirmed_manifest(NULL)
      reviewed_validation_version(-1L)
      invisible(NULL)
    }

    shiny::observeEvent(input$discover, {
      paths <- sarek_parse_path_input(input$paths)
      if (!length(paths)) {
        confirmation_state(NULL)
        invalidate_confirmation()
        status_state("ERROR: Enter at least one absolute server path.")
        validation_state("No manifest has been validated yet.")
        return(invisible(NULL))
      }
      tryCatch({
        table <- sarek_build_discovery_table(
          paths,
          recursive = isTRUE(input$recursive),
          max_files = max_files,
          allowed_roots = allowed_input_roots
        )
        confirmation_state(table)
        sample_review <- sarek_sample_review_table(table)
        selected_sample_state(sarek_sample_key(sample_review$patient_id[[1]], sample_review$sample_id[[1]]))
        invalidate_confirmation()
        recommended_mode <- sarek_recommend_analysis_mode(table)
        if (nzchar(recommended_mode)) {
          shiny::updateSelectInput(session, "analysis_mode", selected = recommended_mode)
        }
        ignored <- as.integer(sarek_shiny_value(attr(table, "ignored_count"), 0L))
        status_state(paste0(
          "ACTION REQUIRED: Discovered ", NROW(table), " supported file", if (NROW(table) == 1L) "" else "s",
          if (ignored > 0L) paste0("; ignored ", ignored, " unsupported file", if (ignored == 1L) "" else "s") else "",
          ". Review every inferred sample field and FASTQ pairing before confirmation",
          if (nzchar(recommended_mode)) paste0(". Suggested analysis mode: ", sarek_analysis_mode_label(recommended_mode)) else "",
          "."
        ))
        validation_state("Discovery complete. The current draft has not been confirmed.")
      }, error = function(error) {
        confirmation_state(NULL)
        invalidate_confirmation()
        status_state(paste0("ERROR: ", conditionMessage(error)))
        validation_state("No manifest has been validated yet.")
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    output$status <- shiny::renderUI({
      message <- status_state()
      kind <- sarek_manifest_status_kind(message)
      shiny::div(
        class = paste("sarek-status-banner", paste0("sarek-status-banner-", kind)),
        shiny::tags$span(message)
      )
    })
    output$summary <- shiny::renderTable(
      sarek_confirmation_summary(confirmation_state()),
      striped = TRUE,
      bordered = TRUE,
      spacing = "s"
    )

    shiny::observeEvent(input$sample_to_edit, {
      selected_sample_state(as.character(sarek_shiny_value(input$sample_to_edit)))
    }, ignoreInit = TRUE)

    selected_sample_key <- shiny::reactive({
      table <- confirmation_state()
      if (is.null(table) || !NROW(table)) return("")
      keys <- unique(sarek_sample_key(table$patient_id, table$sample_id))
      selected <- selected_sample_state()
      if (selected %in% keys) selected else keys[[1]]
    })

    output$readiness <- shiny::renderUI({
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      review <- sarek_sample_review_table(table)
      pairing_issues <- review[startsWith(review$FASTQ_pairing, "Needs attention:"), , drop = FALSE]
      recommended <- sarek_recommend_analysis_mode(table)
      selected_mode <- as.character(sarek_shiny_value(input$analysis_mode))

      pairing_item <- if (!NROW(pairing_issues)) {
        shiny::div(class = "sarek-review-item sarek-review-item-pass", "✓ FASTQ pairing is complete for every included FASTQ sample.")
      } else {
        details <- paste0(pairing_issues$sample_id, " (", sub("^Needs attention: ", "", pairing_issues$FASTQ_pairing), ")")
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          paste0("! FASTQ files need attention: ", paste(details, collapse = "; "), ".")
        )
      }

      mode_item <- if (nzchar(recommended) && identical(selected_mode, recommended)) {
        shiny::div(
          class = "sarek-review-item sarek-review-item-pass",
          paste0("✓ Analysis mode matches the detected roles: ", sarek_analysis_mode_label(recommended), ".")
        )
      } else if (nzchar(recommended)) {
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          paste0(
            "! Detected roles suggest ", sarek_analysis_mode_label(recommended),
            ", but the selected mode is ", sarek_analysis_mode_label(selected_mode), "."
          )
        )
      } else {
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          "! Roles do not support one clear analysis mode yet; review each sample role."
        )
      }

      shiny::div(
        class = "sarek-review-checklist",
        shiny::tags$strong("Pre-validation checklist"),
        pairing_item,
        mode_item
      )
    })

    output$sample_editor <- shiny::renderUI({
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      review <- sarek_sample_review_table(table)
      keys <- sarek_sample_key(review$patient_id, review$sample_id)
      labels <- paste0(review$patient_id, " / ", review$sample_id, " — ", review$role)
      choices <- stats::setNames(keys, labels)
      selected <- selected_sample_key()
      row <- review[keys == selected, , drop = FALSE]
      if (!NROW(row)) row <- review[1, , drop = FALSE]

      normal_ids <- unique(review$sample_id[
        review$patient_id == row$patient_id[[1]] & review$role == "normal"
      ])
      current_match <- as.character(row$matched_normal_id[[1]])
      normal_ids <- unique(c(normal_ids, current_match[nzchar(current_match)]))
      match_choices <- c("No matched normal" = "", stats::setNames(normal_ids, normal_ids))

      shiny::div(
        class = "sarek-sample-editor",
        shiny::h4("Edit one sample"),
        shiny::tags$p(
          class = "muted small-note",
          "Select a sample, correct its metadata once, then apply the change to all associated files."
        ),
        shiny::selectInput(ns("sample_to_edit"), "Sample to review", choices = choices, selected = selected),
        shiny::fluidRow(
          shiny::column(
            6,
            shiny::checkboxInput(ns("edit_include"), "Include this sample", value = isTRUE(row$include[[1]])),
            shiny::textInput(ns("edit_patient_id"), "Patient ID", value = row$patient_id[[1]]),
            shiny::textInput(ns("edit_sample_id"), "Sample ID", value = row$sample_id[[1]])
          ),
          shiny::column(
            6,
            shiny::selectInput(ns("edit_role"), "Sample role", choices = SAREK_SAMPLE_ROLES, selected = row$role[[1]]),
            shiny::selectInput(
              ns("edit_matched_normal_id"),
              "Matched normal (tumor samples only)",
              choices = match_choices,
              selected = current_match
            ),
            shiny::selectInput(
              ns("edit_processing_state"),
              "Processing state",
              choices = SAREK_PROCESSING_STATES,
              selected = row$processing_state[[1]]
            )
          )
        ),
        shiny::tags$p(
          class = "muted small-note",
          paste0("Detected format: ", row$input_format[[1]], ". FASTQ pairing: ", row$FASTQ_pairing[[1]], ".")
        ),
        shiny::actionButton(ns("apply_sample"), "Apply sample changes", class = "btn-primary")
      )
    })

    shiny::observeEvent(input$apply_sample, {
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      patient_id <- trimws(as.character(sarek_shiny_value(input$edit_patient_id)))
      sample_id <- trimws(as.character(sarek_shiny_value(input$edit_sample_id)))
      if (!nzchar(patient_id) || !nzchar(sample_id)) {
        status_state("ERROR: Patient ID and sample ID cannot be blank.")
        return(invisible(NULL))
      }
      tryCatch({
        table <- sarek_apply_sample_update(
          table = table,
          sample_key = selected_sample_key(),
          include = isTRUE(input$edit_include),
          patient_id = patient_id,
          sample_id = sample_id,
          role = input$edit_role,
          matched_normal_id = input$edit_matched_normal_id,
          processing_state = input$edit_processing_state
        )
        confirmation_state(table)
        selected_sample_state(sarek_sample_key(patient_id, sample_id))
        invalidate_confirmation()
        status_state(paste0(
          "ACTION REQUIRED: Updated ", patient_id, " / ", sample_id,
          " across all associated files. Continue reviewing the checklist before validation."
        ))
        validation_state("The draft changed. Validate it again before downloading.")
      }, error = function(error) {
        status_state(paste0("ERROR: ", conditionMessage(error)))
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    output$file_editor <- shiny::renderUI({
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      keys <- sarek_sample_key(table$patient_id, table$sample_id)
      sample_files <- table[keys == selected_sample_key(), , drop = FALSE]
      shiny::req(NROW(sample_files))
      fastq_files <- sample_files[sample_files$input_format == "fastq", , drop = FALSE]
      if (!NROW(fastq_files)) {
        return(shiny::tags$p(class = "muted small-note", "This sample does not use FASTQ input, so lane/read correction is not needed."))
      }

      choices <- stats::setNames(fastq_files$path, basename(fastq_files$path))
      selected_path <- as.character(sarek_shiny_value(input$file_to_edit, fastq_files$path[[1]]))
      if (!selected_path %in% fastq_files$path) selected_path <- fastq_files$path[[1]]
      row <- fastq_files[fastq_files$path == selected_path, , drop = FALSE]
      selected_read <- if (is.na(row$read[[1]])) "" else as.character(row$read[[1]])

      shiny::div(
        class = "sarek-sample-editor",
        shiny::h4("Correct FASTQ lane/read detection"),
        shiny::tags$p(
          class = "muted small-note",
          "Use this only when a valid FASTQ filename was not interpreted correctly. It cannot supply a genuinely missing mate."
        ),
        shiny::selectInput(ns("file_to_edit"), "FASTQ file", choices = choices, selected = selected_path),
        shiny::fluidRow(
          shiny::column(6, shiny::textInput(ns("edit_lane"), "Lane", value = row$lane[[1]], placeholder = "L001")),
          shiny::column(
            6,
            shiny::selectInput(
              ns("edit_read"),
              "Read",
              choices = c("Not detected" = "", "R1" = "1", "R2" = "2"),
              selected = selected_read
            )
          )
        ),
        shiny::actionButton(ns("apply_file_pairing"), "Apply lane/read correction")
      )
    })

    shiny::observeEvent(input$apply_file_pairing, {
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      tryCatch({
        read <- if (nzchar(as.character(sarek_shiny_value(input$edit_read)))) input$edit_read else NA_integer_
        table <- sarek_apply_file_pairing_update(
          table,
          path = input$file_to_edit,
          lane = input$edit_lane,
          read = read
        )
        confirmation_state(table)
        invalidate_confirmation()
        status_state("ACTION REQUIRED: FASTQ lane/read correction applied. Recheck pairing before validation.")
        validation_state("The draft changed. Validate it again before downloading.")
      }, error = function(error) {
        status_state(paste0("ERROR: ", conditionMessage(error)))
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    if (requireNamespace("DT", quietly = TRUE)) {
      output$sample_review_table <- DT::renderDT({
        table <- sarek_sample_review_table(confirmation_state())
        shiny::validate(shiny::need(NROW(table), "Discover inputs to populate this table."))
        table$include <- ifelse(table$include, "☑", "☐")
        names(table) <- c(
          "Include", "Patient ID", "Sample ID", "Role", "Matched normal",
          "Input format", "Processing state", "Files", "FASTQ pairing"
        )
        DT::datatable(
          table,
          rownames = FALSE,
          class = "stripe hover compact nowrap",
          width = "100%",
          selection = list(mode = "single", target = "row"),
          options = list(
            scrollX = TRUE,
            scrollCollapse = FALSE,
            pageLength = 25,
            lengthMenu = list(c(10, 25, 50, 100), c("10", "25", "50", "100")),
            searching = TRUE,
            autoWidth = TRUE
          )
        )
      }, server = FALSE)

      shiny::observeEvent(input$sample_review_table_rows_selected, {
        selected_rows <- input$sample_review_table_rows_selected
        if (!length(selected_rows)) return(invisible(NULL))
        selected_row <- suppressWarnings(as.integer(selected_rows[[1]]))
        review <- sarek_sample_review_table(confirmation_state())
        if (!is.na(selected_row) && selected_row >= 1L && selected_row <= NROW(review)) {
          selected_sample_state(sarek_sample_key(
            review$patient_id[[selected_row]],
            review$sample_id[[selected_row]]
          ))
        }
        invisible(NULL)
      }, ignoreInit = TRUE)

      output$file_detail_table <- DT::renderDT({
        table <- confirmation_state()
        shiny::validate(shiny::need(!is.null(table) && NROW(table), "Discover inputs to populate file details."))
        keys <- sarek_sample_key(table$patient_id, table$sample_id)
        table <- table[keys == selected_sample_key(), c(
          "path", "index", "lane", "read", "size_bytes", "role_confidence", "warning"
        ), drop = FALSE]
        DT::datatable(
          table,
          rownames = FALSE,
          class = "stripe hover compact nowrap",
          width = "100%",
          selection = "none",
          options = list(
            scrollX = TRUE,
            scrollCollapse = FALSE,
            paging = FALSE,
            searching = FALSE,
            info = FALSE,
            autoWidth = TRUE,
            columnDefs = sarek_manifest_column_defs(table)
          )
        )
      }, server = FALSE)
    } else {
      output$sample_review_table <- shiny::renderTable({
        sarek_sample_review_table(confirmation_state())
      }, striped = TRUE, bordered = TRUE, spacing = "xs")
      output$file_detail_table <- shiny::renderTable({
        table <- confirmation_state()
        shiny::req(!is.null(table) && NROW(table))
        keys <- sarek_sample_key(table$patient_id, table$sample_id)
        table[keys == selected_sample_key(), c(
          "path", "index", "lane", "read", "size_bytes", "role_confidence", "warning"
        ), drop = FALSE]
      }, striped = TRUE, bordered = TRUE, spacing = "xs")
    }

    shiny::observeEvent(
      list(
        input$manifest_id,
        input$assay_type,
        input$analysis_mode,
        input$preset,
        input$results_root,
        input$work_root,
        input$species,
        input$assembly,
        input$sarek_genome
      ),
      {
        if (!is.null(confirmed_manifest())) {
          invalidate_confirmation()
          validation_state("The draft settings changed. Validate them again before downloading.")
        }
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(input$confirm, {
      table <- confirmation_state()
      if (is.null(table) || !NROW(table)) {
        confirmed_manifest(NULL)
        validation_state("ERRORS\n- Discover and review inputs before confirming a manifest.")
        status_state("ERROR: Manifest validation could not start. Discover and review inputs first.")
        return(invisible(NULL))
      }
      validation <- sarek_validate_confirmation_table(table)
      mode_errors <- sarek_validate_analysis_mode(table, input$analysis_mode)
      validation_state(sarek_manifest_validation_text(validation, mode_errors))
      if (!isTRUE(validation$valid) || length(mode_errors)) {
        confirmed_manifest(NULL)
        status_state("ERROR: Manifest validation failed. Review the red validation panel and correct every listed item.")
        return(invisible(NULL))
      }
      tryCatch({
        manifest <- sarek_build_manifest(
          confirmation_table = table,
          manifest_id = input$manifest_id,
          created_by = created_by,
          assay_type = input$assay_type,
          analysis_mode = input$analysis_mode,
          preset = input$preset,
          results_root = input$results_root,
          work_root = input$work_root,
          species = input$species,
          assembly = input$assembly,
          sarek_genome = input$sarek_genome,
          allowed_results_roots = allowed_results_roots,
          allowed_work_roots = allowed_work_roots
        )
        confirmed_manifest(manifest)
        validation_version(validation_version() + 1L)
        reviewed_validation_version(-1L)
        validation_state(paste0(
          sarek_manifest_validation_text(validation, mode_errors),
          "\n\nCONFIRMED\n- The JSON manifest is ready to download."
        ))
        if (length(validation$warnings)) {
          status_state("ACTION REQUIRED: Validation passed with warnings. Read the amber validation panel before downloading the JSON manifest.")
        } else {
          status_state("Manifest validation passed. Review the green validation panel before downloading the JSON manifest.")
        }
      }, error = function(error) {
        confirmed_manifest(NULL)
        validation_state(paste0("ERRORS\n- ", conditionMessage(error)))
        status_state("ERROR: Manifest validation failed. Review the red validation panel and correct every listed item.")
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    output$validation <- shiny::renderUI({
      message <- validation_state()
      kind <- sarek_manifest_validation_kind(message)
      shiny::div(
        class = paste("sarek-validation-panel", paste0("sarek-validation-panel-", kind)),
        shiny::div(
          shiny::h4(sarek_manifest_validation_heading(kind)),
          shiny::tags$pre(message)
        )
      )
    })
    output$manifest_preview <- shiny::renderText({
      manifest <- confirmed_manifest()
      if (is.null(manifest)) return("No confirmed manifest is available.")
      text <- jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
      if (nchar(text, type = "bytes") > 30000L) {
        paste0(substr(text, 1L, 30000L), "\n\n[Preview truncated; download contains the complete JSON.]")
      } else {
        text
      }
    })

    output$download_ui <- shiny::renderUI({
      if (is.null(confirmed_manifest())) {
        return(shiny::tags$p(class = "muted small-note", "Confirm a valid draft to enable JSON download."))
      }
      reviewed <- identical(reviewed_validation_version(), validation_version())
      shiny::tagList(
        shiny::checkboxInput(
          session$ns("validation_reviewed"),
          "I reviewed the validation results above",
          value = reviewed
        ),
        if (reviewed) {
          shiny::downloadButton(session$ns("download_manifest"), "Download confirmed JSON")
        } else {
          shiny::tags$p(class = "muted small-note", "Review and acknowledge the validation panel to enable download.")
        }
      )
    })

    shiny::observeEvent(input$validation_reviewed, {
      if (isTRUE(input$validation_reviewed) && !is.null(confirmed_manifest())) {
        reviewed_validation_version(validation_version())
      } else {
        reviewed_validation_version(-1L)
      }
    }, ignoreInit = TRUE)

    output$download_manifest <- shiny::downloadHandler(
      filename = function() {
        paste0(sarek_identifier(input$manifest_id, "sarek_analysis"), ".manifest.json")
      },
      content = function(file) {
        manifest <- confirmed_manifest()
        if (is.null(manifest)) stop("The manifest must be confirmed before download.")
        if (!identical(reviewed_validation_version(), validation_version())) {
          stop("Review and acknowledge the validation results before downloading.")
        }
        sarek_write_manifest(manifest, file)
      },
      contentType = "application/json"
    )

    list(
      confirmation_table = shiny::reactive(confirmation_state()),
      manifest = shiny::reactive(confirmed_manifest()),
      valid = shiny::reactive(!is.null(confirmed_manifest()))
    )
  })
}
