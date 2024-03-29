import pandas as pd
import os
import numpy as np
from abc import ABC, abstractmethod
from scipy.ndimage import uniform_filter1d
import itertools
from scipy.stats import binom


# smoothing truth raw frequency using 1dimension filter
def smooth_freq(df):

    raw_freq = df["truth_freq"].values
    df["smoothed_freq"] = uniform_filter1d(raw_freq, size=7, mode="nearest")
    return df


def load_data(filepath):
    # Check if file exists and continue if not
    if not os.path.exists(filepath):
        return None  
    # read models and add to dict

    raw_pred = pd.read_csv(filepath)
    # print(model,location,date)
    raw_pred["pred_freq"] = raw_pred["median_freq_forecast"]
    if "median_freq_nowcast" in raw_pred.columns:
        raw_pred["pred_freq"] = raw_pred["pred_freq"].fillna(
            value=raw_pred["median_freq_nowcast"]
        )

    if "freq_forecast_lower_50" in raw_pred.columns and "freq_forecast_upper_95" in raw_pred.columns:
        ci_low = raw_pred["freq_forecast_lower_50"]
        ci_high = raw_pred["freq_forecast_upper_95"]
    else:
        ci_low = None
        ci_high = None
    return raw_pred, ci_low, ci_high





# reading truth_set
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

    new_truth = df.merge(truth_set, how="left").fillna(0)

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


# comment
def prep_freq_data(final_set):

    # convert frequencies to arrays
    raw_freq = np.squeeze(final_set[["truth_freq"]].to_numpy())


    if len(raw_freq) == 0:
        return (None,) * 5
    raw_freq[np.isnan(raw_freq)] = 0
    # smoothed freq
    smoothed_freq = np.squeeze(final_set[["smoothed_freq"]].to_numpy())
    smoothed_freq[np.isnan(smoothed_freq)] = 0


    # return seq_total
    seq_count = np.squeeze(final_set[["sequences"]].to_numpy())

    # return seq
    total_seq = np.squeeze(final_set[["total_seq"]].to_numpy())

    # convert predicted frequencies to arrays
    pred_freq = np.squeeze(final_set[["pred_freq"]].to_numpy())

    return raw_freq, pred_freq, seq_count, total_seq, smoothed_freq


def merge_truth_pred(df, location_truth):

    merged_set = pd.merge(df, location_truth, how="left")
    merged_set["sequences"] = merged_set["sequences"].fillna(0)
    # sum sequences of each location and date
    merged_set["total_seq"] = merged_set.groupby(["date", "location"])[
        "sequences"
    ].transform("sum")
    # compute truth frequencies for each variant
    merged_set["truth_freq"] = (
        merged_set["sequences"] / merged_set["total_seq"])
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

    # Calculating error
    # MSE
    mae = MAE()
    # error_df['MAE'] = mae.evaluate(truth_mv_avg,pred_freq)
    error_df["MAE"] = mae.evaluate(smoothed_freq, pred_freq)

    # RMSE
    rmse = RMSE()
    error_df["RMSE"] = rmse.evaluate(smoothed_freq, pred_freq)

    # Logloss error
    logloss = LogLoss()
    error_df["loglik"] = logloss.evaluate(seq_count, total_seq, pred_freq)


    # Adding confidence intervals
    # Computing Coverage
    coverage = Coverage()
    if ci_low is not None and ci_high is not None:
        error_df["freq_forecast_lower_50"] = ci_low
        error_df["freq_forecast_upper_95"] = ci_high

        # Compute coverage
    error_df["coverage"] = Coverage.compute_coverage(smoothed_freq, ci_low, ci_high)
        



    # adding frequencies columns for comparison and diagnostics
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
        #abs_error = prediction - truth
        return abs_error

class Coverage(Scores):
    def __init__(self):
        pass

    def compute_coverage(self, truth, ci_low, ci_high):
        within_interval = (truth >= ci_low) & (truth <= ci_high)
        coverage = within_interval  #Maybe add .mean() here?
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
        # mlr log loss error

    def evaluate(self, seq_value, tot_seq, pred_values):
        loglik = binom.logpmf(k=seq_value, n=tot_seq, p=pred_values)
        return loglik


if __name__ == "__main__":
    locations = [
        "USA",
        "United Kingdom",
        "Brazil",
        "Australia",
        "South Africa",
        "Japan",
        "Malaysia"
    ]
    models = ["GARW", "MLR", "FGA", "Piantham", "dummy"]
    dates = ['2022-01-01', '2022-01-15','2022-02-01','2022-02-15','2022-03-01','2022-03-15',
         '2022-04-01','2022-04-15','2022-05-01','2022-05-15','2022-06-01','2022-06-15',
         '2022-07-01','2022-07-15','2022-08-01','2022-08-15','2022-09-01','2022-09-15',
         '2022-10-01','2022-10-15','2022-11-01','2022-11-15','2022-12-01','2022-12-15']

    # truth_seq_count per variant
    truth_set = load_truthset(
        "../data/time_stamped/truth/seq_counts_truth.tsv"
    )

    # full model output set dict

    score_df_list = []
    # loop thorough different files of model versions
    for model in models:

        for location in locations:
            pred_dic = {}
            # filtering final_truth dataset to run location
            location_truth = truth_set[truth_set["location"] == location]
            # for all specified pivot dates
            for pivot_date in dates:

                filepath = f"../model_estimates_fullRank/cast_estimates_full_{model}/{location}/freq_full_{pivot_date}.csv"

                # Load data
                raw_pred = load_data(filepath)
                print(raw_pred)
                if raw_pred is None:
                    continue

                # Merge predictions and truth set
                merged = merge_truth_pred(raw_pred, location_truth)

                # Make dataframe containing the errors
                error_df = calculate_errors(merged, pivot_date)
                print(error_df)
                if error_df is None:
                    continue

                score_df_list.append(error_df)

    # print(score_df_list)
    score_df = pd.concat(score_df_list)

    # save score output to a csv file
    score_df.to_csv(f"../estimates_fullrank/model_scores_output.csv", index=False)
