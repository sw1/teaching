#!/bin/bash

Rscript -e 'bookdown::render_book("index.Rmd","bookdown::gitbook")'
rm map.tsv* file*.txt file*.html file*.log

exit 0
