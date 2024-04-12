import pandas as pd
import os
import numpy as np
from abc import ABC, abstractmethod
import itertools
from scipy.stats import binom
import yaml


# smoothing truth raw frequency using 1dimension filter
def smooth_freq(df):
    raw_freq = pd.Series(df["truth_freq"])
    smoothed_freq = (
        raw_freq.rolling(window=7, min_periods=1, center=True).mean().values
    )
    df["smoothed_freq"] = smoothed_freq
    return df


def load_data(filepath, sep):
    # Check if file exists and continue if not
    if not os.path.exists(filepath):
        return None  # Return None if file does not exist

    # Read the file into a DataFrame
    raw_pred = pd.read_csv(filepath, sep=sep)

    # Assign values from "median_freq_forecast" column to "pred_freq" column
    raw_pred["pred_freq"] = raw_pred["median_freq_forecast"]

    # Fill missing values in "pred_freq" column with values from "median_freq_nowcast" column if available
    if "median_freq_nowcast" in raw_pred.columns:
        raw_pred["pred_freq"].fillna(
            value=raw_pred["median_freq_nowcast"], inplace=True
        )

    # Fixing column name if using old runs
    if "freq_nowcast_upper_95" in raw_pred.columns:
        raw_pred["freq_upper_95"] = raw_pred["freq_nowcast_upper_95"]

    # Check if forecasting confidence intervals are available in the DataFrame
    if (
        "freq_forecast_lower_95" in raw_pred.columns
        and "freq_forecast_upper_95" in raw_pred.columns
    ):
        # Add columns for confidence intervals to the DataFrame
        raw_pred["ci_low"] = raw_pred["freq_forecast_lower_95"]
        raw_pred["ci_high"] = raw_pred["freq_forecast_upper_95"]
        if "freq_lower_95" in raw_pred.columns:
            raw_pred["ci_low"] = raw_pred[["ci_low", "freq_lower_95"]].min(
                axis=1
            )
        if "freq_upper_95" in raw_pred.columns:
            raw_pred["ci_high"] = raw_pred[
                ["ci_high", "freq_upper_95"]
            ].max(axis=1)

    # Return the modified DataFrame raw_pred
    return raw_pred


# Loading retrospective frequencies
def load_truthset(path):
    truth_set = pd.read_csv(path, sep="\t")

    all_dates = pd.unique(truth_set["date"])
    all_loc = pd.unique(truth_set["location"])
    all_var = pd.unique(truth_set["variant"])

    combined = [all_dates, all_loc, all_var]
    df = pd.DataFrame(
        columns=["date", "location", "variant"],
        data=list(itertools.product(*combined)),
    )

    new_truth = df.merge(truth_set, how="left")#.fillna(0)

    new_truth["total_seq"] = new_truth.groupby(["date", "location"])[
        "sequences"
    ].transform("sum")

    new_truth["truth_freq"] = new_truth["sequences"] / new_truth["total_seq"]
    new_truth = new_truth.sort_values(by=["location", "variant", "date"])
    new_truth = (
        new_truth.groupby(["location", "variant"], group_keys=False)
        .apply(smooth_freq)
        .reset_index(drop=True)
    )
    return new_truth


def prep_freq_data(final_set):

    # Convert frequencies to arrays
    raw_freq = np.squeeze(final_set[["truth_freq"]].to_numpy(), axis=-1)

    if len(raw_freq) == 0:
        return (None,) * 7
    raw_freq[np.isnan(raw_freq)] = 0

    # Convert smoothed frequencies to arrays
    smoothed_freq = np.squeeze(final_set[["smoothed_freq"]].to_numpy())
    smoothed_freq[np.isnan(smoothed_freq)] = 0

    # Convert sequences and total sequnces to arrays
    seq_count = np.squeeze(final_set[["sequences"]].to_numpy())
    total_seq = np.squeeze(final_set[["total_seq"]].to_numpy())

    # Convert predicted frequencies to arrays
    pred_freq = np.squeeze(final_set[["pred_freq"]].to_numpy())

    # Convert credible intervals to arrays
    if "ci_low" in final_set.columns:
        ci_low = np.squeeze(final_set[["ci_low"]].to_numpy())
    else:
        ci_low = None

    if "ci_high" in final_set.columns:
        ci_high = np.squeeze(final_set[["ci_high"]].to_numpy())
    else:
        ci_high = None

    return (
        raw_freq,
        pred_freq,
        seq_count,
        total_seq,
        smoothed_freq,
        ci_low,
        ci_high,
    )


def merge_truth_pred(df, location_truth):

    merged_set = pd.merge(df, location_truth, how="left")
    merged_set["sequences"] = merged_set["sequences"]#.fillna(0)

    # Compute total sequences for each location and date
    merged_set["total_seq"] = merged_set.groupby(["date", "location"])[
        "sequences"
    ].transform("sum")

    # Compute retrospective frequencies for each variant
    merged_set["truth_freq"] = (
        merged_set["sequences"] / merged_set["total_seq"]
    )
    return merged_set[merged_set["pred_freq"].notnull()]


def calculate_errors(merged, pivot_date):
    # Compute model dates and leads from pivot_date
    model_dates = pd.to_datetime(merged["date"])
    lead = (model_dates - pd.to_datetime(pivot_date)).dt.days

    error_df = pd.DataFrame(
        {
            "location": location,
            "model": model,
            "pivot_date": pd.to_datetime(pivot_date),
            "lead": lead,
            "variant": merged["variant"],
        }
    )

    # unpacking prepped_data values
    (
        raw_freq,
        pred_freq,
        seq_count,
        total_seq,
        smoothed_freq,
        ci_low,
        ci_high,
    ) = prep_freq_data(merged)

    if raw_freq is None:
        return None

    # Computing metrics
    # MSE
    mae = MAE()
    error_df["MAE"] = mae.evaluate(smoothed_freq, pred_freq)

    # RMSE
    rmse = RMSE()
    error_df["RMSE"] = rmse.evaluate(smoothed_freq, pred_freq)

    # Logloss error
    logloss = LogLoss()
    error_df["loglik"] = logloss.evaluate(seq_count, total_seq, pred_freq)

    # Computing Coverage
    coverage = Coverage()
    if ci_low is not None and ci_high is not None:
        error_df["coverage"] = coverage.compute_coverage(
            smoothed_freq, ci_low, ci_high
        )

    # Adding frequencies columns for comparison and diagnostics
    error_df["total_seq"] = total_seq
    error_df["raw_freq"] = raw_freq
    error_df["smoothed_freq"] = smoothed_freq
    error_df["pred_freq"] = pred_freq
    error_df["date"] = model_dates
    return error_df


class Scores(ABC):
    @abstractmethod
    def __init__(self):
        pass


class MAE(Scores):
    def __init__(self):
        pass

    def evaluate(self, truth, prediction):
        abs_error = np.abs(truth - prediction)
        return abs_error


class Coverage(Scores):
    def __init__(self):
        pass

    def compute_coverage(self, truth, ci_low, ci_high):
        within_interval = (truth >= ci_low) & (truth <= ci_high)
        coverage = within_interval.astype(int)
        return coverage


class MSE(Scores):
    def __init__(self):
        pass

    def evaluate(self, truth, prediction):
        squared_error = np.square(truth - prediction)
        return squared_error


class RMSE(Scores):
    def __init__(self):
        pass

    def evaluate(self, truth, prediction):
        squared_error = np.square(truth - prediction)
        return np.sqrt(squared_error)


class LogLoss(Scores):
    def __init__(self):
        pass

    def evaluate(self, seq_value, tot_seq, pred_values):
        loglik = binom.logpmf(k=seq_value, n=tot_seq, p=pred_values)
        return loglik


if __name__ == "__main__":
    with open("../config.yaml", "r") as config:
        config = yaml.safe_load(config)

    dates = config["main"]["estimation_dates"]
    locations = config["main"]["locations"]
    models = config["main"]["models"]

    # Retrospective sequence counts
    truth_set = load_truthset(
        "../data/time_stamped/truth/seq_counts_truth.tsv"
    )

    # Loading model predictions and computing errors
    score_df_list = []
    for model in models:
        for location in locations:
            pred_dic = {}
            # Filtering to location of increast
            location_truth = truth_set[truth_set["location"] == location]
            if location_truth is None:
                print(
                    f"{location_truth} is not in the retrospective frequencies."
                )
                continue

            for pivot_date in dates:
                filepath = f"../estimates/{model}/{location}/freq_full_{pivot_date}.tsv"
                # if no .tsv, check for .csv
                if not os.path.exists(filepath):
                    filepath = f"../estimates/{model}/{location}/freq_full_{pivot_date}.csv"

                # Load data
                sep = "\t" if filepath.endswith(".tsv") else ","
                raw_pred = load_data(filepath, sep=sep)
                if raw_pred is None:
                    print(
                        f"No {model} predictions found for location {location} on analysis date {pivot_date}"
                    )
                    continue

                # Merge predictions and truth set
                merged = merge_truth_pred(raw_pred, location_truth)

                # Make dataframe containing the errors
                error_df = calculate_errors(merged, pivot_date)
                if error_df is None:
                    continue

                score_df_list.append(error_df)
    score_df = pd.concat(score_df_list)

    # Save score output to a tsv file
    error_filepath = "../errors/"
    if not os.path.exists(error_filepath):
        os.makedirs(error_filepath)
    score_df.to_csv("../errors/model_scores.tsv", sep="\t", index=False)
