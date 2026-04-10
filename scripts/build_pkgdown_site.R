#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgdown)
})

pkg <- "."

pkgdown::build_home(pkg = pkg, quiet = FALSE)
pkgdown::build_reference(pkg = pkg, lazy = FALSE, preview = FALSE, devel = FALSE)
pkgdown::build_news(pkg = pkg)

# pkgdown::build_articles() always uses the callr-backed article path in the
# current pkgdown release. Building articles directly keeps the same output
# format while avoiding that wrapper and its broken error handling.
article_section <- pkgdown:::section_init(pkg, "articles")
pkgdown:::build_articles_index(article_section)

article_names <- article_section$vignettes$name[article_section$vignettes$type == "rmd"]
for (name in article_names) {
  pkgdown:::build_article(
    name,
    pkg = pkg,
    lazy = FALSE,
    seed = 1014L,
    new_process = FALSE,
    quiet = FALSE
  )
}

pkgdown::build_search(pkg = pkg)
