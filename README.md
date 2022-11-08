# ncov-forecasting-fit

**Eslam Abousamra** <sup>1,2</sup>, **Marlin Figgins** <sup>1,3,4</sup>, **Trevor Bedford** <sup>1,4</sup>

<sup>1</sup> *Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA* <br />
<sup>2</sup> *Department of Epidemiology, University of Washington, Seattle, WA, USA* 
<sup>3</sup> *Department of Applied Mathematics, University of Washington, Seattle, WA, USA*
<sup>4</sup> *Howard Hughes Medical Institute, Seattle, WA, USA*




The **#ncov-forecasting-fit** repository hosts a data curation and a live-forecasting framework to process pathogen variant data at time-stamped intervals. The framework is built to standardize and estimate accurate real-time nowcast and forecast **targets** and to facilitate comparisons of forecasting and nowcasting accuracy between different statistical models. Using the framework, the purpose of the study is to work with live surveillance data to investigate the empirical side of evolutionary forecasting including growth advantages, frequencies estimates, and cases and to provide a scoring framework of different modelling approaches.



### Abstract










### Data processing
Outline pipeline going from sequences to estimates 



### Model comparison and evaluation (Analysis I)

Data: format data in the form of sequence counts per location per model per observation data (as known of that date)
what specific time period and frequency to run these analysis on?
Modeling: Naive, Piantham, MLR, FGA, GARW



### Scoring Script

MSE, MAE, logLoss
Computation of truth set vs estimates. How to estimate the truth set, should be logloss, mlr run or other?
Visualization: Violin plots, scatterplots, (?) of errors over time. Group error to show variants errors over time independently.
Conclusions: Inferences and discussion of findings 


### Comparison of growth advantages (Analysis II)
Data: format data in the form of growth advantages estimates per location per model per observation date.
Explore growth advantages change over time relative to the first analysis
Visualization: discuss how to we want to visualize the variation of growth advantage per variant
Conclusion: inferences and discussion of findings

### Examination of data quality reflection to error (Analysis III)
Data: Submission delays, missing data (sequence counts per date), errors
Visualization: Scatterplot of sequence counts per date and errors
Analysis: infer relationship between sequence efforts and quality to model errors (Fit model)
Conclusion: Make statements regarding the relationship






