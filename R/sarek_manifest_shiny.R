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
    shiny::fluidRow(
      shiny::column(
        4,
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
        shiny::textAreaInput(
          ns("paths"),
          "Input files or folders",
          value = "",
          rows = 8,
          placeholder = "/absolute/path/sample_T_L001_R1.fastq.gz\n/absolute/path/sample_T_L001_R2.fastq.gz"
        ),
        shiny::checkboxInput(ns("recursive"), "Search inside subfolders", value = FALSE),
        shiny::actionButton(ns("discover"), "Discover inputs", class = "btn-primary"),
        shiny::br(),
        shiny::br(),
        shiny::verbatimTextOutput(ns("status"))
      ),
      shiny::column(
        8,
        shiny::tags$div(
          class = "read-source-note",
          shiny::tags$strong("Step 2: corroborate the guesses"),
          shiny::tags$p(
            "Edit include, patient/sample IDs, role, matched normal, format, processing state, lane, or read."
          ),
          shiny::tags$p(
            class = "muted small-note",
            "File paths, detected indexes, sizes, confidence, and warnings are read-only."
          )
        ),
        shiny::tableOutput(ns("summary")),
        sarek_manifest_table_output(ns("confirmation_table"))
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
        shiny::textInput(ns("results_root"), "Results root", value = ""),
        shiny::textInput(ns("work_root"), "Work root", value = ""),
        shiny::tags$p(
          class = "muted small-note",
          "Both must be absolute server paths. This step records them but does not create them."
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
    shiny::verbatimTextOutput(ns("validation")),
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
    confirmation_state <- shiny::reactiveVal(NULL)
    confirmed_manifest <- shiny::reactiveVal(NULL)
    status_state <- shiny::reactiveVal(
      "Enter one or more server paths, then select Discover inputs."
    )
    validation_state <- shiny::reactiveVal(
      "No manifest has been validated yet."
    )

    shiny::observe({
      shiny::updateTextInput(session, "results_root", value = default_results_root)
      shiny::updateTextInput(session, "work_root", value = default_work_root)
    })

    invalidate_confirmation <- function() {
      confirmed_manifest(NULL)
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
        invalidate_confirmation()
        ignored <- as.integer(sarek_shiny_value(attr(table, "ignored_count"), 0L))
        status_state(paste0(
          "Discovered ", NROW(table), " supported file", if (NROW(table) == 1L) "" else "s",
          if (ignored > 0L) paste0("; ignored ", ignored, " unsupported file", if (ignored == 1L) "" else "s") else "",
          ". Review the table before confirmation."
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

    output$status <- shiny::renderText(status_state())
    output$summary <- shiny::renderTable(
      sarek_confirmation_summary(confirmation_state()),
      striped = TRUE,
      bordered = TRUE,
      spacing = "s"
    )

    if (requireNamespace("DT", quietly = TRUE)) {
      output$confirmation_table <- DT::renderDT({
        table <- confirmation_state()
        shiny::validate(shiny::need(!is.null(table) && NROW(table), "Discover inputs to populate this table."))
        read_only <- setdiff(names(table), sarek_editable_columns())
        disabled <- match(read_only, names(table)) - 1L
        DT::datatable(
          table,
          rownames = FALSE,
          editable = list(target = "cell", disable = list(columns = disabled)),
          options = list(
            scrollX = TRUE,
            pageLength = 25,
            lengthMenu = list(c(10, 25, 50, 100), c("10", "25", "50", "100")),
            autoWidth = FALSE
          )
        )
      }, server = FALSE)

      shiny::observeEvent(input$confirmation_table_cell_edit, {
        table <- sarek_apply_confirmation_edit(
          confirmation_state(),
          input$confirmation_table_cell_edit
        )
        confirmation_state(table)
        invalidate_confirmation()
        validation_state("The draft changed. Validate it again before downloading.")
      }, ignoreInit = TRUE)
    } else {
      output$confirmation_table <- shiny::renderTable({
        table <- confirmation_state()
        shiny::validate(shiny::need(!is.null(table) && NROW(table), "Discover inputs to populate this table."))
        table
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
        return(invisible(NULL))
      }
      validation <- sarek_validate_confirmation_table(table)
      mode_errors <- sarek_validate_analysis_mode(table, input$analysis_mode)
      validation_state(sarek_manifest_validation_text(validation, mode_errors))
      if (!isTRUE(validation$valid) || length(mode_errors)) {
        confirmed_manifest(NULL)
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
        validation_state(paste0(
          sarek_manifest_validation_text(validation, mode_errors),
          "\n\nCONFIRMED\n- The JSON manifest is ready to download."
        ))
      }, error = function(error) {
        confirmed_manifest(NULL)
        validation_state(paste0("ERRORS\n- ", conditionMessage(error)))
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    output$validation <- shiny::renderText(validation_state())
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
      shiny::downloadButton(session$ns("download_manifest"), "Download confirmed JSON")
    })

    output$download_manifest <- shiny::downloadHandler(
      filename = function() {
        paste0(sarek_identifier(input$manifest_id, "sarek_analysis"), ".manifest.json")
      },
      content = function(file) {
        manifest <- confirmed_manifest()
        if (is.null(manifest)) stop("The manifest must be confirmed before download.")
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
