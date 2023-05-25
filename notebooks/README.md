



## Running the Models

This notebook provides code for variant frequency forecasting and growth advantage estimation using different models ['GARW','MLR', 'Piantham', 'FGA'] and provides specification
for a naive model. It does that by performing the following functions:

1. Specification of analysis period and geographical settings: Define the observation dates and geographical settings for the analysis.

2. Assigning Parameters and forecasting period: Set the parameters and forecasting periods for the models along with model definition specification.

4. Running the Models to obtain variant frequency forecasts and growth advantage estimates for each model and location.

This process can be found in the notebook ./models_run.ipynb.


## Creating data sets as of observation date
This generates the sequence counts for variants of interest from available sequence data as of a given date as well as an additional retrospective truth data set i.e. including all sequences up to now.

This process can be found in the notebook ./creating-data-sets.ipynb.
