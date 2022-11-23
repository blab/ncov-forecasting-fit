# Manuscript

Rebuild with `Rake`. This will run:

```
pdflatex -draftmode ncov-forecasting-fit
pdflatex -draftmode ncov-forecasting-fit
bibtex ncov-forecasting-fit
pdflatex ncov-forecasting-fit
```

to generate the PDF `ncov-forecasting-fit.pdf`. However, on subsequent builds it will skip steps if they are not required.
The PDF `ncov-forecasting-fit.pdf` is intentionally not versioned.

