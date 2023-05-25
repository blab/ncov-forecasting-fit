# Scripts

## Computing model errors based on truth sets

This script computes model errors for different evolutionary forecasting models ("GARW", "MLR", "FGA", "Piantham", "dummy") based on provided truth sets by performing following functions:

1. Loads the truth set data.
2. Processes the model output files for each location and model.
3. Merges the model predictions with the corresponding truth set.
4. Calculates error metrics such as Mean Absolute Error (MAE), Root Mean Squared Error (RMSE), and Log Loss.
5. Generates an output DataFrame containing the error metrics and additional frequency columns.
6. Saves the error metrics to a CSV file.

## Formatting growth advantages from `evofr` models
