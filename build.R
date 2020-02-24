#!/usr/bin/env Rscript

bookdown::render_book(".", "bookdown::gitbook", new_session=TRUE)
bookdown::render_book(".", "bookdown::pdf_book", new_session=TRUE)
