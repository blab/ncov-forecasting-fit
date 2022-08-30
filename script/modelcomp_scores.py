from statistics import mode
import pandas as pd
import os
import numpy as np
from abc import ABC, abstractmethod
from scipy.ndimage import uniform_filter1d
import itertools


#smoothing truth raw frequency using 1dimension filter
def smooth_freq(df):

    raw_freq = df['truth_freq']
    df['smoothed_freq'] = uniform_filter1d(raw_freq, size=7, mode = "nearest")

    return df


#reading truth_set
def load_truthset(path):

    truth_set = pd.read_csv(path, sep="\t")

    all_dates = pd.unique(truth_set['date'])
    all_loc = pd.unique(truth_set['location'])
    all_var = pd.unique(truth_set['variant'])
    
    combined = [all_dates, all_loc, all_var]
    df = pd.DataFrame(columns= ['date', 'location', 'variant'], 
        data = list(itertools.product(*combined)))

    new_truth = df.merge(truth_set, how = "left").fillna(0)

    new_truth['total_seq'] = new_truth.groupby(['date', 'location'])['sequences'].transform('sum')

    new_truth['truth_freq'] = new_truth['sequences']/new_truth['total_seq']
    new_truth = new_truth.sort_values(by=["location", "variant", "date"])
    new_truth = new_truth.groupby(['location', 'variant']).apply(smooth_freq).reset_index(drop=True)
    
    return new_truth



def prep_freq_data(final_set):
    
    #convert frequencies to arrays
    truth_values = np.squeeze(final_set[['smoothed_freq']].to_numpy())
    if len(truth_values) == 0: 
        return (None,)*5

    raw_freq = np.squeeze(final_set[['truth_freq']].to_numpy())
    #return seq_total 
    seq_value = np.squeeze(final_set[['sequences']].to_numpy())
    #return seq
    total_seq =np.squeeze(final_set[['total_seq']].to_numpy())

    #convert predicted frequencies to arrays
    pred_values =np.squeeze(final_set[['pred_freq']].to_numpy())


    return raw_freq, pred_values, seq_value, total_seq, truth_values

def merge_truth_pred(df, location_truth):
    merged_set = pd.merge(df, location_truth, how = 'left')
    return merged_set[merged_set["pred_freq"].notnull()]




class Scores(ABC):
    @abstractmethod
    def __init__(self):
        pass
class MAE(Scores):
    def __init__(self):
        pass
    def evaluate(self,truth, prediction):
        abs_error = np.abs(truth - prediction)
        return abs_error
        
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

class LogLoss(Scores): #to-do
    def __init__(self):
        pass
        #mlr log loss error
    def evaluate(self, truth, prediction):
        pass


if __name__=='__main__':
    locations = ["USA","Japan", "United Kingdom"]
    models = ["GARW", "MLR", "FGA"]
    dates = ['2022-04-15','2022-04-22','2022-04-29','2022-05-06',
         '2022-05-13','2022-05-20','2022-05-27','2022-06-03',
         '2022-06-10','2022-06-17','2022-06-24','2022-06-30']

    #reading truth_sequence count dataset per variant
    truth_set = load_truthset("../data/truth/seq_counts_truth.tsv")



    #full model output set dict
    final_sets = {}

    #loop thorough different files of model versions
    for model in models:

        for location in locations:
            pred_dic = {}
            #filtering final_truth dataset to run location
            location_truth = truth_set[truth_set['location']==location]
            for date in dates:
                
                filepath = f"../plot-est/cast_estimates_full_{model}/{location}/freq_full_{date}.csv"
                
                #Check if file exists and continue if not
                if not os.path.exists(filepath):
                    continue

                #read models and add to dict
                raw_pred = pd.read_csv(filepath)
                #adding predicted frequency column
                raw_pred['pred_freq'] = raw_pred['median_freq_forecast'].fillna(raw_pred['median_freq_nowcast'])
                

                pred_dic[date] =  raw_pred
            
            #merging predicted frequency with smoothed frequency
            final_sets_location = {k: merge_truth_pred(df,location_truth) for k,df in pred_dic.items()} 


            final_sets[location, model] = final_sets_location

        
            #final_sets_new  =  final_sets[location, model].fillna('NA')

            
            #final_sets_new =  {k: fillna(v) for k, v in final_sets[location, model].values()}



    score_df_list = []

    for model in models:

        for location in locations:
            #prep_Freq_data through conversion to arrays to pass in the metric class
            prepped_data = {k: prep_freq_data(v) for k,v in final_sets[location, model].items()}

            for k, v in prepped_data.items():
                model_dates = pd.to_datetime(final_sets[location, model][k]['date'])
                pivot_date = pd.to_datetime(k)
                lead = (model_dates-pivot_date).dt.days
                variants = final_sets[location, model][k]['variant']

                #creating error dict with results
                error_df = pd.DataFrame({'location': location, 'model': model, 'pivot_date': pivot_date, 'lead': lead, 'variant':variants})

                #unpacking prepped_data values
                raw_freq, pred_values, seq_value, total_seq, truth_values = v 
            
                
                if truth_values is None:
                    continue

                #Calculating MAE error
                mae = MAE()  
                error_df['MAE'] = mae.evaluate(truth_values,pred_values)


                #Calculating RMSE error
                rmse = RMSE()
                error_df['RMSE'] = rmse.evaluate(truth_values,pred_values)

                #adding frequencies columns for comparison and diagnostics
                error_df['raw_freq'] = raw_freq
                error_df['smoothed_freq'] = truth_values
                error_df['pred_freq'] = pred_values
                error_df['date'] = model_dates


                
                score_df_list.append(error_df)



                
    score_df = pd.concat(score_df_list)


    #save score output to a csv file
    score_df.to_csv(f"../estimates/model_scores_output3.csv",index = False)


