#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgdown)
})

pkg <- "."
expected_site_default <-
  'document.body.classList.add("preset-interaction")'
site_script <- paste(
  readLines(file.path(pkg, "pkgdown", "extra.js"), warn = FALSE),
  collapse = "\n"
)
if (!grepl(expected_site_default, site_script, fixed = TRUE)) {
  stop(
    "pkgdown/extra.js does not activate the interaction preset by default.",
    call. = FALSE
  )
}

pkgdown::build_home(pkg = pkg, quiet = FALSE)
pkgdown::build_reference(pkg = pkg, lazy = FALSE, preview = FALSE, devel = FALSE)
pkgdown::build_news(pkg = pkg)

# pkgdown::build_articles() always uses the callr-backed article path in the
# current pkgdown release. Building articles directly keeps the same output
# format while avoiding that wrapper and its broken error handling.
article_section <- pkgdown:::section_init(pkg, "articles")
pkgdown:::build_articles_index(article_section)

article_names <- article_section$vignettes$name[article_section$vignettes$type == "rmd"]
expected_theme_classes <-
  'document.body.classList.add("palette-red","preset-interaction")'

for (name in article_names) {
  pkgdown:::build_article(
    name,
    pkg = pkg,
    lazy = FALSE,
    seed = 1014L,
    new_process = FALSE,
    quiet = FALSE
  )

  article_path <- file.path(pkg, "docs", "articles", paste0(name, ".html"))
  article_html <- paste(readLines(article_path, warn = FALSE), collapse = "\n")
  if (!grepl(expected_theme_classes, article_html, fixed = TRUE)) {
    stop(
      "Article '", name,
      "' did not activate the expected red interaction Albers theme.",
      call. = FALSE
    )
  }
}

pkgdown::build_search(pkg = pkg)
