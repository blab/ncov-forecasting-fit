<<<<<<< HEAD
# Compiling

The LaTeX manuscript can be compiled by running `rake` from this directory. This will run:
```
pdflatex -draftmode ...
bibtex ...
pdflatex -draftmode ...
pdflatex ...
```
or skip steps when possible.
=======
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

>>>>>>> 9e83a793b72db391a6823d211e00a08aed8f21c8
