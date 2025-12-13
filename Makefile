.PHONY: coverage test check vignette docs

coverage:
	@echo "Running coverage via covr..."
	@Rscript -e "if (!requireNamespace('covr', quietly = TRUE)) { stop('Package \'covr\' is required. Install it with install.packages(\'covr\').') }" \
		-e "cov <- covr::package_coverage(type = c('tests')); print(cov); covr::report(cov, file = 'coverage.html'); cat('\nCoverage report written to coverage.html\n')"

test:
	@echo "Running testthat suite..."
	@Rscript -e "testthat::test_dir('tests/testthat', reporter = 'summary')"

check:
	@echo "Running R CMD check..."
	@R CMD build . >/dev/null 2>&1; \
	TARBALL=$$(ls -1t *.tar.gz | head -n1); \
	R CMD check "$$TARBALL" --no-manual --as-cran

vignette:
	@echo "Rendering vignette(s)..."
	@Rscript -e "rmarkdown::render('vignettes/intro-gdsfmri.Rmd', quiet = TRUE)"

docs: vignette
