# ncov-forecasting-fit

**Eslam Abousamra** <sup>1,2</sup>, **Marlin Figgins** <sup>1,3</sup>, **Trevor Bedford** <sup>1,4</sup>

<sup>1</sup> *Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA*
<sup>2</sup> *Department of Epidemiology, University of Washington, Seattle, WA, USA* 
<sup>3</sup> *Department of Applied Mathematics, University of Washington, Seattle, WA, USA*
<sup>4</sup> *Howard Hughes Medical Institute, Seattle, WA, USA*




The **#ncov-forecasting-fit** repository hosts a data curation and a live-forecasting framework to process pathogen variant data at time-stamped intervals. The framework is built to standardize and estimate accurate real-time nowcast and forecast **targets** and to facilitate comparisons of forecasting and nowcasting accuracy between different statistical models. Using the framework, the purpose of the study is to work with live surveillance data to investigate the empirical side of evolutionary forecasting including growth advantages, frequencies estimates, and cases and to provide a scoring framework of different modelling approaches.



### Abstract










### Data

Script to obtain raw data can be found `script/get_versioned_data.py`, please refer to `script/README.md` for more information on how to obtain the data.


### Running the models

`notebooks/auto_renewal_accuracy.ipynb` notebook contains code to generate time-stamped estimates of variant frequencies, growth advantages, and cases for a specified number of countries using five different models that vary in complexity. 



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


