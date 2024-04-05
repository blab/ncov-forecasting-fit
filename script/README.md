# Scripts

## Computing model errors based on truth sets

This script computes model errors for different evolutionary forecasting models ("GARW", "MLR", "FGA", "Piantham", "dummy") based on provided retrospective sequence counts.

The script does this in the following steps:
1. Loads the retrospective sequence counts to produce a truth set of retrospective smooth frequencies.
2. Loads the model predictions from each model and estimation date.
3. Merges the model predictions with the corresponding truth set.
4. Calculates error metrics such as Mean Absolute Error (MAE), Root Mean Squared Error (RMSE), and coverage.
5. Generates an output DataFrame containing the error metrics and additional retrospective frequency columns.
6. Exports the error metrics to a file `../estimates/model_scores.tsv`.

This process can be found in the script `./compute_model_scores.py`.

## Formatting growth advantages from `evofr` models

This script formats growth advantages from evolutionary forecasting models for different locations and models. It calculates the growth advantage of variants for all specified pivot dates.
It does that by performing the following functions:

1. Formats the growth advantage data by creating a DataFrame with columns for location, model, pivot date, and variant.
2. Defines the list of locations and models to process.
3. Calls the formatting function to create a formatted DataFrame for each combination of model, location, and pivot date.
4. Exports the formatted growth advantage data for all locations, models, and variants for each pivot date to a file `../estimates/ga_output.csv`.

This process can be found in the script `./tidy_growth_adv.py`.
