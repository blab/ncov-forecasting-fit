import pandas as pd
import os

locations = [
    "USA",
    "Japan",
    "United Kingdom",
    "Brazil",
    "Australia",
    "South Africa",
]
models = ["GARW", "MLR", "FGA", "Piantham"]
dates = [
    "2022-01-01",
    "2022-01-15",
    "2022-02-01",
    "2022-02-15",
    "2022-03-01",
    "2022-03-15",
    "2022-04-01",
    "2022-04-15",
    "2022-05-01",
    "2022-05-15",
    "2022-06-01",
    "2022-06-15",
    "2022-07-01",
    "2022-07-15",
    "2022-08-01",
    "2022-08-15",
    "2022-09-01",
    "2022-09-15",
    "2022-10-01",
    "2022-10-15",
    "2022-11-01",
    "2022-11-15",
    "2022-12-01",
    "2022-12-15",
]


def format_growth_advantages(ga_data, pivot_date, model, location):
    ga_df = pd.DataFrame(
        {
            "location": location,
            "model": model,
            "pivot_date": pd.to_datetime(pivot_date),
            "variant": ga_data["variant"],
        }
    )

    # unpacking prepped_data values
    if model == "GARW":
        ga_df["date"] = ga_data["date"]

    value_columns = ["median_ga", "ga_upper_80", "ga_lower_80"]
    ga_df[value_columns] = ga_data[value_columns]
    return ga_df


def main(dates, models, locations):

    # full model output set dict
    ga_formatted_list = []

    # loop thorough different files of model versions
    for model in models:

        for location in locations:
            # filtering final_truth dataset to run location
            for pivot_date in dates:

                filepath = (
                    f"../estimates/{model}/{location}"
                    + f"/full_growth_advantages_{pivot_date}.csv"
                )

                # Check if file exists and continue if not
                if not os.path.exists(filepath):
                    continue

                # read models and add to dict
                ga_data = pd.read_csv(filepath)
                ga_formatted = format_growth_advantages(
                    ga_data, pivot_date, model, location
                )

                # creating ga dict
                ga_formatted_list.append(ga_formatted)

    # CONCAT ALL DATASETS
    ga_formatted = pd.concat(ga_formatted_list)

    # save score output to a csv file
    # output csv of all datasets for each country
    ga_formatted.to_csv("../estimates/ga_output.csv", index=False)
    return None


if __name__ == "__main__":
    main(dates, models, locations)
