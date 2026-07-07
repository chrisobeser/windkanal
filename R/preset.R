#' Kalibrierte Setting-Presets
#'
#' Ein Preset ist eine **Behauptung ueber die Realitaet**: ein
#' benannter Parametersatz fuer [sim_stream()], bei dem jede Zahl
#' eine Quelle traegt (Attribut `sources`; Status: `belegt` oder
#' `offen`). Unkalibrierte Groessen bleiben bewusst draussen --
#' die Stichproben-Groessen (`n_therapists`, `patients_per_therapist`,
#' `weeks_accrual`) sind Skalen-Entscheidungen des Nutzers, keine
#' Realitaets-Behauptungen.
#'
#' Verfuegbar:
#' * `"ambulanz_de"` -- deutsche Hochschul-/Ausbildungsambulanz:
#'   - `n_sessions = 24`: Kurzzeittherapie-Kontingent; 77% der
#'     ambulanten Richtlinien-Therapien sind KZT (<= 24 Sitzungen),
#'     Hochschulambulanzen bieten typisch 24 Sitzungen an
#'     (Quellen: Aerzteblatt-Versorgungsstudie; HH-Hochschulambulanz;
#'     Abruf 2026-07-03).
#'   - `icc = 0.05`: Therapeuten-Anteil an der Outcome-Varianz in
#'     **Routineversorgung**. Wichtig: In kontrollierten Studien ist
#'     der ICC nahe null (Baldwin et al., 2011, Cogn Behav Ther,
#'     40(1): Median -.026 ueber 495 Schaetzungen -- verifiziert);
#'     die ~5-9%-Werte stammen aus naturalistischen Settings
#'     (Kanon-Zitation dafuer noch zu verifizieren).
#'   - `dropout = 0.017`: pro Sitzung, so gewaehlt, dass kumulativ
#'     ~33% abbrechen -- die KODAP-Abbruchquote
#'     (Quelle: methoden/Margraf_2021_KODAP_Cooperative_Revolution.pdf).
#'   - `dropout_informative = 0`: **offen** -- Staerke der Kopplung
#'     Abbruch~Befinden ist noch unkalibriert.
#'   - `shape = "loglinear"`, `mean_slope = -0.40`: Besserungskurve
#'     negativ beschleunigt (Dose-Response-Forschung, Howard-
#'     Tradition; log-lineare Modelle schlagen lineare konsistent).
#'     Kalibrierung: Prae-Post-Effektstaerke in Routineversorgung
#'     ~d=0.9 (Meta-Analyse routinemaessig erbrachter Therapien:
#'     Depression d=.96, Angst d=.80; Abruf 2026-07-03). Rechnung:
#'     d * SD_quer (~1.43) / log(24) = 0.9*1.43/3.18 ~= 0.40.
#'   - `sd_slope = 0.25`: Verlaufs-Heterogenitaet -- **offen**
#'     (plausibel gesetzt, unkalibriert).
#'   - `reliability_score = 0.90`: uebliche interne Konsistenz
#'     gaengiger Symptommasse (~alpha .9) -- Kanon-Zitat offen.
#'
#' Weitere Presets:
#' * `"allianz_beierl2021"` -- der Allianz-Prozess nach Beierl et al.
#'   (2021, *Frontiers in Psychiatry*, 12:602648; N = 230,
#'   PTSD-Kognitive-Therapie -- **Anker, nicht Kanon**): Traegheit
#'   `alliance_ar = 0.75`, Kopplung Allianz->Symptome
#'   `coupling = -0.13` (therapeuten-geratet), Gegenkopplung
#'   `coupling_reverse = -0.12`. Kombinierbar mit dem Setting-Preset:
#'   `c(preset("ambulanz_de"), preset("allianz_beierl2021"))`.
#'
#' @param name Preset-Name: `"ambulanz_de"` oder
#'   `"allianz_beierl2021"`.
#' @return Benannte Liste von [sim_stream()]-Argumenten mit Attribut
#'   `sources` (data.frame: Parameter, Wert, Quelle, Status).
#'
#' @examples
#' p <- preset("ambulanz_de")
#' attr(p, "sources")
#' s <- do.call(sim_stream, c(p, list(n_therapists = 10,
#'                                    patients_per_therapist = 8,
#'                                    seed = 1)))
#'
#' @export
preset <- function(name = "ambulanz_de") {
  if (name == "allianz_beierl2021") {
    args <- list(
      alliance = TRUE,
      alliance_ar = 0.75,
      coupling = -0.13,
      coupling_reverse = -0.12
    )
    attr(args, "sources") <- data.frame(
      parameter = c("alliance_ar", "coupling", "coupling_reverse"),
      wert      = c("0.75", "-0.13", "-0.12"),
      quelle    = rep(paste("Beierl et al. (2021), Frontiers in",
                            "Psychiatry 12:602648; N=230, PTSD --",
                            "Anker, nicht Kanon"), 3),
      status    = c("belegt (AR .75-.79)",
                    "belegt (therapeuten-geratet; patient. n.s.)",
                    "belegt"),
      stringsAsFactors = FALSE
    )
    return(args)
  }
  if (name != "ambulanz_de") {
    stop("Unbekanntes Preset: '", name,
         "'. Verfuegbar: 'ambulanz_de', 'allianz_beierl2021'.",
         call. = FALSE)
  }
  args <- list(
    n_sessions = 24,
    icc = 0.05,
    dropout = 0.017,
    dropout_informative = 0,
    shape = "loglinear",
    mean_slope = -0.40,
    sd_slope = 0.25,
    icc_slope = 0.17,
    reliability_score = 0.90
  )
  attr(args, "sources") <- data.frame(
    parameter = c("n_sessions", "icc", "dropout",
                  "dropout_informative", "shape", "mean_slope",
                  "sd_slope", "icc_slope", "reliability_score"),
    wert      = c("24", "0.05", "0.017/Sitzung (~33% kumulativ)", "0",
                  "loglinear", "-0.40 (Ziel: prae-post d~0.9)",
                  "0.25", "0.17", "0.90"),
    quelle    = c(
      "KZT-Kontingent; 77% KZT (Aerzteblatt-Studie); HH-Ambulanz typ. 24",
      "Johns et al. 2019 (ClinPsychRev): naturalistisch Oe 5% (0.2-21); Trials strittig (Baldwin 2011: ~0 vs Johns: 8.2%)",
      "KODAP ~33% Abbrueche (Margraf et al. 2021, PDF in quellen/)",
      "keine belastbare Zahl gefunden",
      "Dose-Response: log-linear schlaegt linear (Howard-Tradition)",
      "Gaskell et al. 2022, Adm Policy Ment Health 50(1): d=.96/.80",
      "plausibel gesetzt",
      "Lutz et al. 2007: ~17% der Raten-Varianz auf Therapeuten-Ebene",
      "uebliche interne Konsistenz gaengiger Symptommasse"
    ),
    status    = c("belegt", "belegt (naturalistisch, Johns et al. 2019)",
                  "belegt", "offen", "belegt",
                  "belegt (Kalibrier-Rechnung dokumentiert)",
                  "offen",
                  "belegt (Anker; Journal-Zitat zu verifizieren)",
                  "plausibel (Kanon-Zitat offen)"),
    stringsAsFactors = FALSE
  )
  args
}
