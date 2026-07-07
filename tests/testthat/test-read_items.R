test_that("read_items extrahiert Struktur, aber keine Texte", {
  p <- testthat::test_path("fixtures", "mini_battery.csv")
  it <- read_items(p)
  expect_equal(nrow(it), 4)           # note + heading uebersprungen
  expect_named(it, c("item", "scale", "n_options", "reversed", "rev_marker"))
  expect_false(any(grepl("niedergeschlagen|verstanden",
                         unlist(it))))  # kein Itemtext im Output!
  expect_equal(it$scale, c("sym", "sym", "all", "sym"))  # R-Suffix gestrippt
  expect_equal(it$n_options, c(4L, 4L, 2L, 4L))
  expect_equal(it$rev_marker, c(FALSE, TRUE, FALSE, FALSE))
  expect_equal(it$reversed, c(FALSE, FALSE, FALSE, TRUE))  # formr-R-Suffix
  expect_equal(attr(it, "source_file"), "mini_battery.csv")
})

test_that("Nicht-formr-Dateien werden klar abgelehnt", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(a = 1), tmp, row.names = FALSE)
  expect_error(read_items(tmp), "formr")
})
