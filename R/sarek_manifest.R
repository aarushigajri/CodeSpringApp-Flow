# Pure helpers for discovering inputs and building a user-confirmed Sarek manifest.
# This module deliberately has no Shiny or pipeline-runtime dependency.

SAREK_MANIFEST_SCHEMA_VERSION <- "1.0"
SAREK_INSPECTOR_VERSION <- "0.1.0"
SAREK_SUPPORTED_INPUT_FORMATS <- c("fastq", "ubam", "bam", "cram", "vcf", "bcf")
SAREK_SAMPLE_ROLES <- c("germline", "tumor", "normal", "unknown")
SAREK_PROCESSING_STATES <- c(
  "unknown",
  "unmapped",
  "aligned",
  "duplicate_marked",
  "recalibrated",
  "analysis_ready",
  "variant_calls"
)

sarek_text <- function(value, default = "") {
  if (is.null(value) || !length(value) || is.na(value[[1]])) return(default)
  value <- trimws(as.character(value[[1]]))
  if (nzchar(value)) value else default
}

sarek_identifier <- function(value, fallback = "sample") {
  value <- gsub("[^A-Za-z0-9_.-]+", "_", sarek_text(value, fallback))
  value <- gsub("^[_.-]+|[_.-]+$", "", value)
  if (nzchar(value)) substr(value, 1L, 128L) else fallback
}

sarek_detect_input_format <- function(path) {
  name <- tolower(basename(sarek_text(path)))
  if (grepl("\\.(fastq|fq)(\\.gz)?$", name, perl = TRUE)) return("fastq")
  if (grepl("\\.ubam$", name, perl = TRUE)) return("ubam")
  if (grepl("\\.bam$", name, perl = TRUE)) return("bam")
  if (grepl("\\.cram$", name, perl = TRUE)) return("cram")
  if (grepl("\\.vcf(\\.gz)?$", name, perl = TRUE)) return("vcf")
  if (grepl("\\.bcf$", name, perl = TRUE)) return("bcf")
  ""
}

sarek_normalize_existing_path <- function(path) {
  path <- sarek_text(path)
  if (!nzchar(path)) stop("An input path is empty.")
  if (!file.exists(path) && !dir.exists(path)) stop("Input path does not exist: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sarek_find_index <- function(path, format = sarek_detect_input_format(path)) {
  candidates <- switch(
    format,
    bam = c(paste0(path, ".bai"), sub("\\.bam$", ".bai", path, ignore.case = TRUE)),
    cram = c(paste0(path, ".crai"), sub("\\.cram$", ".crai", path, ignore.case = TRUE)),
    vcf = c(paste0(path, ".tbi"), paste0(path, ".csi")),
    bcf = c(paste0(path, ".csi")),
    character(0)
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates)) normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE) else ""
}

sarek_discover_input_paths <- function(paths, recursive = FALSE, max_files = 5000L) {
  paths <- unique(vapply(paths, sarek_normalize_existing_path, character(1)))
  max_files <- suppressWarnings(as.integer(max_files))
  if (!is.finite(max_files) || max_files < 1L) stop("max_files must be a positive integer.")

  candidates <- character(0)
  for (path in paths) {
    if (dir.exists(path)) {
      candidates <- c(
        candidates,
        list.files(
          path,
          full.names = TRUE,
          recursive = isTRUE(recursive),
          all.files = FALSE,
          include.dirs = FALSE,
          no.. = TRUE
        )
      )
    } else {
      candidates <- c(candidates, path)
    }
  }

  candidates <- unique(candidates[file.exists(candidates) & !dir.exists(candidates)])
  formats <- vapply(candidates, sarek_detect_input_format, character(1))
  supported <- nzchar(formats)
  supported_paths <- candidates[supported]

  if (length(supported_paths) > max_files) {
    stop(
      "Input discovery found ", length(supported_paths),
      " supported files, above the safety limit of ", max_files,
      ". Select a narrower folder or increase the administrator limit."
    )
  }

  list(
    source_paths = paths,
    supported_paths = supported_paths,
    ignored_paths = candidates[!supported]
  )
}

sarek_filename_tokens <- function(path, format = sarek_detect_input_format(path)) {
  name <- basename(path)
  stem <- name
  stem <- sub("\\.(fastq|fq)(\\.gz)?$", "", stem, ignore.case = TRUE, perl = TRUE)
  stem <- sub("\\.(ubam|bam|cram|vcf|bcf)(\\.gz)?$", "", stem, ignore.case = TRUE, perl = TRUE)

  lane <- if (grepl("(^|[_.-])L[0-9]{3}([_.-]|$)", stem, ignore.case = TRUE, perl = TRUE)) {
    sub(".*(?:^|[_.-])L([0-9]{3})(?:[_.-]|$).*", "L\\1", stem, ignore.case = TRUE, perl = TRUE)
  } else {
    ""
  }

  read <- if (grepl("(^|[_.-])R1([_.-]|$)", stem, ignore.case = TRUE, perl = TRUE)) {
    1L
  } else if (grepl("(^|[_.-])R2([_.-]|$)", stem, ignore.case = TRUE, perl = TRUE)) {
    2L
  } else {
    NA_integer_
  }

  sample_stem <- stem
  sample_stem <- sub("([_.-])R[12]([_.-][0-9]{3})?$", "", sample_stem, ignore.case = TRUE, perl = TRUE)
  sample_stem <- sub("([_.-])L[0-9]{3}$", "", sample_stem, ignore.case = TRUE, perl = TRUE)
  sample_stem <- sub("([_.-])S[0-9]+$", "", sample_stem, ignore.case = TRUE, perl = TRUE)
  sample_stem <- sarek_identifier(sample_stem)

  role <- "unknown"
  role_confidence <- "none"
  if (grepl("(^|[_.-])tumou?r([_.-]|$)", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "tumor"
    role_confidence <- "medium"
  } else if (grepl("(^|[_.-])normal([_.-]|$)", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "normal"
    role_confidence <- "medium"
  } else if (grepl("(^|[_.-])germline([_.-]|$)", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "germline"
    role_confidence <- "medium"
  } else if (grepl("([_.-])T$", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "tumor"
    role_confidence <- "low"
  } else if (grepl("([_.-])N$", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "normal"
    role_confidence <- "low"
  }

  patient_stem <- sample_stem
  patient_stem <- gsub(
    "(^|[_.-])(tumou?r|normal|germline)([_.-]|$)",
    "_",
    patient_stem,
    ignore.case = TRUE,
    perl = TRUE
  )
  if (role %in% c("tumor", "normal") && identical(role_confidence, "low")) {
    patient_stem <- sub("([_.-])[TN]$", "", patient_stem, ignore.case = TRUE, perl = TRUE)
  }
  patient_stem <- sarek_identifier(patient_stem, sample_stem)

  processing_state <- switch(
    format,
    fastq = "unmapped",
    ubam = "unmapped",
    vcf = "variant_calls",
    bcf = "variant_calls",
    "unknown"
  )

  warning <- if (identical(role, "unknown")) {
    "Sample role requires confirmation."
  } else {
    paste0("Sample role was inferred from the filename with ", role_confidence, " confidence.")
  }

  list(
    sample_id = sample_stem,
    patient_id = patient_stem,
    role = role,
    role_confidence = role_confidence,
    lane = lane,
    read = read,
    processing_state = processing_state,
    warning = warning
  )
}

sarek_build_discovery_table <- function(paths, recursive = FALSE, max_files = 5000L) {
  discovery <- sarek_discover_input_paths(paths, recursive = recursive, max_files = max_files)
  if (!length(discovery$supported_paths)) {
    stop("No supported FASTQ, uBAM, BAM, CRAM, VCF, or BCF files were found.")
  }

  rows <- lapply(discovery$supported_paths, function(path) {
    format <- sarek_detect_input_format(path)
    tokens <- sarek_filename_tokens(path, format)
    info <- file.info(path)
    data.frame(
      include = TRUE,
      patient_id = tokens$patient_id,
      sample_id = tokens$sample_id,
      role = tokens$role,
      input_format = format,
      processing_state = tokens$processing_state,
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      index = sarek_find_index(path, format),
      lane = tokens$lane,
      read = tokens$read,
      size_bytes = as.numeric(info$size[[1]]),
      role_confidence = tokens$role_confidence,
      warning = tokens$warning,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  attr(table, "source_paths") <- discovery$source_paths
  attr(table, "ignored_count") <- length(discovery$ignored_paths)
  table
}

sarek_validate_confirmation_table <- function(table) {
  errors <- character(0)
  warnings <- character(0)
  required <- c(
    "include",
    "patient_id",
    "sample_id",
    "role",
    "input_format",
    "processing_state",
    "path"
  )

  missing <- setdiff(required, names(table))
  if (length(missing)) {
    return(list(
      valid = FALSE,
      errors = paste("Missing required columns:", paste(missing, collapse = ", ")),
      warnings = warnings
    ))
  }

  selected <- table[as.logical(table$include), , drop = FALSE]
  if (!NROW(selected)) errors <- c(errors, "At least one input file must be included.")
  if (!NROW(selected)) return(list(valid = FALSE, errors = errors, warnings = warnings))

  identifier_pattern <- "^[A-Za-z0-9][A-Za-z0-9_.-]*$"
  bad_patient <- !grepl(identifier_pattern, selected$patient_id, perl = TRUE)
  bad_sample <- !grepl(identifier_pattern, selected$sample_id, perl = TRUE)
  if (any(bad_patient)) errors <- c(errors, "Patient IDs may contain only letters, numbers, periods, underscores, and hyphens.")
  if (any(bad_sample)) errors <- c(errors, "Sample IDs may contain only letters, numbers, periods, underscores, and hyphens.")
  if (any(!selected$role %in% SAREK_SAMPLE_ROLES)) errors <- c(errors, "One or more sample roles are invalid.")
  if (any(!selected$input_format %in% SAREK_SUPPORTED_INPUT_FORMATS)) errors <- c(errors, "One or more input formats are invalid.")
  if (any(!selected$processing_state %in% SAREK_PROCESSING_STATES)) errors <- c(errors, "One or more processing states are invalid.")
  if (anyDuplicated(selected$path)) errors <- c(errors, "The same input path is included more than once.")
  if (any(!file.exists(selected$path))) errors <- c(errors, "One or more selected input files no longer exist.")
  if (any(file.access(selected$path, mode = 4) != 0)) errors <- c(errors, "One or more selected input files are not readable.")

  sample_key <- paste(selected$patient_id, selected$sample_id, sep = "::")
  formats_by_sample <- split(selected$input_format, sample_key)
  mixed <- names(formats_by_sample)[vapply(formats_by_sample, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(mixed)) errors <- c(errors, "A sample cannot mix different input formats.")

  if (any(selected$role == "unknown")) warnings <- c(warnings, "Every unknown sample role must be confirmed before pipeline submission.")
  missing_index <- selected$input_format %in% c("bam", "cram", "vcf", "bcf") &
    (!"index" %in% names(selected) || !nzchar(selected$index))
  if (any(missing_index)) warnings <- c(warnings, "One or more indexed formats do not currently have a detected index.")

  list(valid = !length(errors), errors = unique(errors), warnings = unique(warnings))
}

sarek_validate_analysis_mode <- function(table, analysis_mode) {
  selected <- table[as.logical(table$include), , drop = FALSE]
  errors <- character(0)
  roles <- unique(selected$role)
  formats <- unique(selected$input_format)

  if (identical(analysis_mode, "germline") && any(!roles %in% "germline")) {
    errors <- c(errors, "Germline mode requires every included sample to have the germline role.")
  }
  if (identical(analysis_mode, "tumor_only") && any(!roles %in% "tumor")) {
    errors <- c(errors, "Tumor-only mode requires every included sample to have the tumor role.")
  }
  if (identical(analysis_mode, "matched_tumor_normal")) {
    by_patient <- split(selected$role, selected$patient_id)
    incomplete <- names(by_patient)[!vapply(
      by_patient,
      function(x) any(x == "tumor") && any(x == "normal"),
      logical(1)
    )]
    if (length(incomplete)) {
      errors <- c(errors, paste0(
        "Matched tumor-normal mode requires at least one tumor and one normal for each patient: ",
        paste(incomplete, collapse = ", ")
      ))
    }
  }
  if (identical(analysis_mode, "annotation_only") && any(!formats %in% c("vcf", "bcf"))) {
    errors <- c(errors, "Annotation-only mode accepts only VCF or BCF inputs.")
  }
  if (!analysis_mode %in% c("germline", "tumor_only", "matched_tumor_normal", "annotation_only")) {
    errors <- c(errors, "The analysis mode is invalid.")
  }

  unique(errors)
}

sarek_manifest_file_record <- function(row) {
  result <- list(path = as.character(row$path[[1]]))
  if ("index" %in% names(row) && nzchar(sarek_text(row$index[[1]]))) result$index <- as.character(row$index[[1]])
  if ("lane" %in% names(row) && nzchar(sarek_text(row$lane[[1]]))) result$lane <- as.character(row$lane[[1]])
  if ("read" %in% names(row) && !is.na(row$read[[1]])) result$read <- as.integer(row$read[[1]])
  result
}

sarek_build_manifest <- function(
  confirmation_table,
  manifest_id,
  created_by,
  assay_type,
  analysis_mode,
  preset,
  results_root,
  work_root,
  source_paths = attr(confirmation_table, "source_paths"),
  species = "human",
  assembly = "GRCh38",
  sarek_genome = "GATK.GRCh38"
) {
  validation <- sarek_validate_confirmation_table(confirmation_table)
  mode_errors <- sarek_validate_analysis_mode(confirmation_table, analysis_mode)
  errors <- unique(c(validation$errors, mode_errors))
  if (length(errors)) stop(paste(errors, collapse = "\n"))

  selected <- confirmation_table[as.logical(confirmation_table$include), , drop = FALSE]
  results_root <- normalizePath(results_root, winslash = "/", mustWork = FALSE)
  work_root <- normalizePath(work_root, winslash = "/", mustWork = FALSE)
  if (!startsWith(results_root, "/") || !startsWith(work_root, "/")) {
    stop("Results and work roots must be absolute server paths.")
  }

  patient_ids <- unique(selected$patient_id)
  patients <- lapply(patient_ids, function(patient_id) {
    patient_rows <- selected[selected$patient_id == patient_id, , drop = FALSE]
    sample_ids <- unique(patient_rows$sample_id)
    samples <- lapply(sample_ids, function(sample_id) {
      sample_rows <- patient_rows[patient_rows$sample_id == sample_id, , drop = FALSE]
      notes <- unique(sample_rows$warning[nzchar(sample_rows$warning)])
      sample <- list(
        sample_id = sample_id,
        role = unique(sample_rows$role)[[1]],
        input_format = unique(sample_rows$input_format)[[1]],
        processing_state = unique(sample_rows$processing_state)[[1]],
        files = lapply(seq_len(NROW(sample_rows)), function(i) sarek_manifest_file_record(sample_rows[i, , drop = FALSE]))
      )
      if (length(notes)) sample$confirmation_notes <- notes
      sample
    })

    relationships <- list()
    if (identical(analysis_mode, "matched_tumor_normal")) {
      tumors <- unique(patient_rows$sample_id[patient_rows$role == "tumor"])
      normals <- unique(patient_rows$sample_id[patient_rows$role == "normal"])
      for (tumor in tumors) {
        for (normal in normals) {
          relationships[[length(relationships) + 1L]] <- list(
            type = "matched_tumor_normal",
            sample_ids = c(tumor, normal)
          )
        }
      }
    }

    list(
      patient_id = patient_id,
      samples = samples,
      relationships = relationships
    )
  })

  source_paths <- unique(vapply(
    if (is.null(source_paths) || !length(source_paths)) selected$path else source_paths,
    function(path) normalizePath(path, winslash = "/", mustWork = FALSE),
    character(1)
  ))

  list(
    schema_version = SAREK_MANIFEST_SCHEMA_VERSION,
    manifest_id = sarek_identifier(manifest_id, "sarek_analysis"),
    status = "confirmed",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    created_by = sarek_text(created_by, "unknown"),
    reference = list(
      species = species,
      assembly = assembly,
      sarek_genome = sarek_genome
    ),
    assay = list(
      type = assay_type,
      intervals = NULL
    ),
    analysis = list(
      mode = analysis_mode,
      preset = sarek_text(preset, "core")
    ),
    storage = list(
      results_root = results_root,
      work_root = work_root
    ),
    patients = patients,
    provenance = list(
      source_paths = source_paths,
      inspector_version = SAREK_INSPECTOR_VERSION
    )
  )
}

sarek_write_manifest <- function(manifest, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite R package is required to write a Sarek manifest.")
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".sarek_manifest_", tmpdir = dirname(path), fileext = ".json")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  jsonlite::write_json(
    manifest,
    path = temporary,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )
  if (!file.rename(temporary, path)) stop("Could not atomically save the Sarek manifest: ", path)
  path
}
