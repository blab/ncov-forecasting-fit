# Fitness models provide accurate short-term forecasts of SARS-CoV-2 variant frequency

**Eslam Abousamra**<sup>1,2</sup>, **Marlin Figgins**<sup>1,3</sup>, **Trevor Bedford**<sup>1,2,4</sup>

<sup>1</sup> *Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA*
<sup>2</sup> *Department of Epidemiology, University of Washington, Seattle, WA, USA*
<sup>3</sup> *Department of Applied Mathematics, University of Washington, Seattle, WA, USA*
<sup>4</sup> *Howard Hughes Medical Institute, Seattle, WA, USA*

The **#ncov-forecasting-fit** repository hosts a data curation and a live-forecasting framework to process pathogen variant data at time-stamped intervals. The framework is built to standardize and estimate accurate real-time nowcast and forecast **targets** and to facilitate comparisons of forecasting and nowcasting accuracy between different statistical models. Using the framework, the purpose of the study is to work with live surveillance data to investigate the empirical side of evolutionary forecasting including growth advantages, frequencies estimates, and cases and to provide a scoring framework of different modelling approaches.

### Abstract

Genomic surveillance of pathogen evolution is essential for public health response, treatment strategies, and vaccine development. In the context of SARS-COV-2, multiple models have been developed including Multinomial Logistic Regression (MLR), Fixed Growth Advantage (FGA), Growth Advantage Random Walk (GARW) and Piantham that use observed variant sequence counts through time to analyze variant dynamics. These models provide estimates of variant fitness and can be used to forecast changes in variant frequency. We introduce a framework for evaluating real-time forecasts of variant frequencies, and apply this framework to the evolution of SARS-CoV-2 during 2022 in which multiple new viral variants emerged and rapidly spread through the population. We compare models across representative countries with different intensities of genomic surveillance. Retrospective assessment of model accuracy highlights that most models of variant frequency perform well and are able to produce reasonable forecasts. We find that the simple MLR model provides less than 1% mean absolute error when forecasting 30 days out for countries with robust genomic surveillance. We investigate impacts of sequence quantity and quality across countries on forecast accuracy and conduct systematically downsampling to identify that 1000 sequences per week is fully sufficient for accurate short-term forecasts. We conclude that fitness models represent a useful prognostic tool for short-term evolutionary forecasting.

### Data

Notebook to obtain raw data can be found `notebooks/creating-data-sets.ipynb`, please refer to `notebooks/README.md` for more information on how to obtain the data.


### Running the models

`notebooks/models_run.ipynb` notebook contains code to generate time-stamped estimates of variant frequencies, growth advantages, and cases for a specified number of countries using five different models that vary in complexity.



### I. Model comparison and evaluation

**Input**: formatted data in the form of sequence counts per location per model per observation data (as known of that date)  <br />
**Models**: Naive, Piantham, MLR, FGA, GARW

### Scoring Script

Scoring estimates for different models can be generated using `script/modelcomp_scores.py`


### II. Comparison of growth advantages
**Input**: format data in the form of growth advantages estimates per location per model per observation date.
GA estimates for different models can be generated using `script/tidy_growth_adv.py`

### III. Examination of data quality reflection to error

Effect of emergence of variants on model estimation.


### Figures

Figures can be generated using `script/script_result_vis.rmd`
