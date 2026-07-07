#' formr-Itemtabelle einlesen -- nur die Struktur
#'
#' Liest eine formr-Fragebogentabelle (.xlsx oder .csv) und extrahiert
#' **ausschliesslich die Struktur**: Itemname, Skalenzuordnung,
#' Antwortformat, Umkehr-Marker. Die Item-*Texte* (`label`) werden
#' bewusst **nicht** uebernommen -- Itemtexte sind oft urheberrechtlich
#' geschuetzt; die Struktur ist es nicht. So koennen Struktur-Dateien
#' bedenkenlos geteilt werden, waehrend Inhalts-Dateien privat bleiben
#' (Lizenz-Disziplin als Code, nicht als Policy).
#'
#' Konventionen (aus den realen formr-Tabellen des Projekts):
#' * Items = Zeilen, deren `type` mit `mc` beginnt (ausser
#'   `mc_heading`); `note`/`submit`-Zeilen werden uebersprungen.
#' * `[rev]`-Marker in `explanations` werden als `rev_marker`
#'   uebernommen. **Achtung, Datei-Konvention:** die Bedeutung ist
#'   projektspezifisch (in den BCE-Tabellen dieses Projekts =
#'   *revidierte Fassung*, nicht Umkehr-Kodierung!).
#' * Skala = Namens-Praefix ohne Item-Nummer (`ctq_3` -> `ctq`,
#'   `BCE_01` -> `bce`); ein `R`-Suffix nach der Nummer wird dabei
#'   ebenfalls entfernt.
#' * **Offizielle formr-Umkehr-Konvention:** Items, deren Name nach
#'   der Nummer auf `R` endet (`bfi_2R`), werden bei der formr-
#'   Skalenaggregation automatisch umgepolt -> Spalte `reversed`.
#' * Antwortoptionen, in dieser Reihenfolge gesucht: eigenes
#'   `choices`-Blatt der Datei (formr-Standard: Spalte `list_name`,
#'   eine Zeile pro Option) -> `choice*`-Spalten der Zeile bzw. der
#'   zugehoerigen `mc_heading`-Zeile -> bekannte Kuerzel
#'   (`ja_nein` = 2) -> sonst `NA` mit Warnung.
#'
#' @param path Pfad zur formr-Tabelle (`.xlsx` oder `.csv`).
#' @return `data.frame` mit einer Zeile pro Item: `item`, `scale`,
#'   `n_options`, `reversed` (offizielle formr-R-Suffix-Konvention),
#'   `rev_marker` (projektspezifische [rev]-Kennzeichnung). Attribut
#'   `source_file` fuer Provenienz.
#' @export
read_items <- function(path) {
  stopifnot(file.exists(path))
  choice_lists <- integer(0)
  raw <- if (grepl("[.]xlsx?$", path, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Fuer .xlsx wird das Paket 'readxl' gebraucht.", call. = FALSE)
    }
    # choices-Blatt (formr-Standard) einsammeln, falls vorhanden
    sheets <- readxl::excel_sheets(path)
    for (sh in sheets[-1]) {
      cand <- as.data.frame(readxl::read_excel(path, sheet = sh,
                                               col_types = "text"))
      names(cand) <- tolower(trimws(names(cand)))
      if ("list_name" %in% names(cand)) {
        # Excel-Konvention: list_name oft nur in der ersten Zeile
        # der Liste -> nach unten auffuellen (LOCF)
        ln <- trimws(cand$list_name)
        ln[!is.na(ln) & ln == ""] <- NA
        for (k in seq_along(ln)) {
          if (is.na(ln[k]) && k > 1) ln[k] <- ln[k - 1]
        }
        tab <- table(ln[!is.na(ln)])
        choice_lists <- c(choice_lists,
                          stats::setNames(as.integer(tab), names(tab)))
      }
    }
    as.data.frame(readxl::read_excel(path, col_types = "text"))
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE,
                    colClasses = "character")
  }
  names(raw) <- tolower(trimws(names(raw)))
  if (!all(c("type", "name") %in% names(raw))) {
    stop("Keine formr-Tabelle: Spalten 'type' und 'name' fehlen.",
         call. = FALSE)
  }
  typ <- trimws(ifelse(is.na(raw$type), "", raw$type))

  choice_cols <- grep("^choice", names(raw), value = TRUE)
  n_choices_in_row <- function(i) {
    if (!length(choice_cols)) return(0L)
    v <- unlist(raw[i, choice_cols])
    sum(!is.na(v) & trimws(v) != "")
  }

  # Antwortlisten aus mc_heading-Zeilen einsammeln (Name -> n Optionen)
  heads <- which(grepl("^mc_heading", typ))
  head_opts <- stats::setNames(
    vapply(heads, n_choices_in_row, integer(1)),
    trimws(raw$name[heads])
  )
  bekannt <- c(ja_nein = 2L)

  items <- which(grepl("^mc", typ) & !grepl("^mc_heading", typ))
  if (!length(items)) {
    stop("Keine Item-Zeilen (type 'mc*') gefunden.", call. = FALSE)
  }

  n_opt <- vapply(items, function(i) {
    eigene <- n_choices_in_row(i)
    if (eigene > 0) return(as.integer(eigene))
    # Referenz auf Antwortliste: "mc ctq_choices" / "mc_button ja_nein"
    ref <- trimws(sub("^mc[a-z_]*", "", typ[i]))
    if (ref %in% names(choice_lists)) return(choice_lists[[ref]])
    if (ref %in% names(head_opts) && head_opts[[ref]] > 0)
      return(head_opts[[ref]])
    if (ref %in% names(bekannt))  return(bekannt[[ref]])
    NA_integer_
  }, integer(1))
  if (anyNA(n_opt)) {
    warning("Antwortoptionen nicht bestimmbar fuer: ",
            paste(raw$name[items][is.na(n_opt)], collapse = ", "))
  }

  expl <- if ("explanations" %in% names(raw)) raw$explanations else ""
  nm <- trimws(raw$name[items])
  out <- data.frame(
    item      = nm,
    scale     = tolower(sub("_?[0-9]+_?R?$", "", nm)),
    n_options = n_opt,
    reversed  = grepl("[0-9]_?R$", nm),
    rev_marker = grepl("\\[rev\\]", ifelse(is.na(expl[items]), "",
                                           expl[items]))
  )
  rownames(out) <- NULL
  attr(out, "source_file") <- basename(path)
  out
}
