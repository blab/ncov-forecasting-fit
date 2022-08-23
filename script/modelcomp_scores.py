from itertools import dropwhile
from pickle import TRUE
from statistics import mode
from typing import final
import pandas as pd
import sklearn
from sklearn.metrics import log_loss, brier_score_loss
import os
import numpy as np
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
from datetime import datetime




#comment
def prep_freq_data(final_set):
    
    #convert frequencies to arrays
    truth_values = np.squeeze(final_set[['truth_freq']].to_numpy())
    if len(truth_values) == 0: 
        return (None,)*5
    #smoothing frequencies
    #smoothing
    filter_length = 7
    truth_mv_avg = np.convolve(truth_values, np.ones((filter_length)), mode = 'same')
    truth_mv_avg /= filter_length
    
    
    #return seq_total 
    seq_value = np.squeeze(final_set[['sequences']].to_numpy())
    #return seq
    total_seq =np.squeeze(final_set[['total_seq']].to_numpy())

    #convert predicted frequencies to arrays
    values =np.squeeze(final_set[['pred_freq']].to_numpy())

    return truth_values, values, seq_value, total_seq, truth_mv_avg



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
    models = ["GARW", "MLR", "FGA", "GARW-N"]
    dates = ['2022-01-24', '2022-02-04','2022-02-08','2022-02-18','2022-02-23',
         '2022-02-28','2022-03-03','2022-03-08','2022-03-15',
         '2022-03-21','2022-03-25','2022-04-07','2022-04-14','2022-04-27'
         ,'2022-05-06','2022-05-17','2022-05-20','2022-05-28','2022-06-09'
         ,'2022-06-14','2022-06-22']
    #Latest model run "truth"

    #truth_seq_count per variant
    truth_set = pd.read_csv("../data/2022-06-30/seq_counts_2022-06-30.tsv", sep="\t")
    #sum sequences of each location and date
    truth_set['total_seq'] = truth_set.groupby(['date', 'location'])['sequences'].transform('sum')
    #compute truth frequencies for each variant
    truth_set['truth_freq'] = truth_set['sequences']/truth_set['total_seq']



    #full model output set dict


    final_sets = {}

    #loop thorough different files of model versions
    for model in models:

        for location in locations:
            pred_dic = {}
            #filtering final_truth dataset to run location
            location_truth = truth_set[truth_set['location']==location]
            location_truth = location_truth[['date','location','variant','truth_freq', 'total_seq','sequences']]
            #print(location_truth)
            for date in dates:
                
                filepath = f"../plot-est/cast_estimates_full_{model}/{location}/freq_full_{date}.csv"
                
                #Check if file exists and continue if not
                if not os.path.exists(filepath):
                    continue
                    #read models and add to dict
                
                raw_pred = pd.read_csv(filepath)
                
                raw_pred['pred_freq'] = raw_pred['median_freq_forecast'].fillna(raw_pred['median_freq_nowcast'])
            
     
                pred_dic[date] =  raw_pred


            final_sets_location = {k: pd.merge(location_truth,df) for k,df in pred_dic.items()}            

            final_sets[location, model] = final_sets_location


    #print(final_sets['USA', 'GARW']['2022-05-17'])

    test = prep_freq_data(final_sets['USA','GARW']['2022-05-17'])

    #print(final_sets[location,model].values())


    score_df_list = []

    for model in models:

        for location in locations:
            #prep_Freq_data for truth values
            prepped_data = {k: prep_freq_data(v) for k,v in final_sets[location, model].items()}

            for k, v in prepped_data.items():
                model_dates = pd.to_datetime(final_sets[location, model][k]['date'])
                pivot_date = pd.to_datetime(k)
                lead = (model_dates-pivot_date).dt.days
                variants = final_sets[location, model][k]['variant']


                error_df = pd.DataFrame({'location': location, 'model': model, 'pivot_date': pivot_date, 'lead': lead, 'variant':variants})
                
                #unpacking prepped_data values
                true_freq, pred_freq, sequences, total_sequences, truth_mv_avg = v 
                #one of sets likely empty causing error
                if true_freq is None: 
                    continue

                #Calculating error
                mae = MAE()  
                error_df['MAE'] = mae.evaluate(truth_mv_avg,pred_freq)
                #print(v[1])
                #MSE
                rmse = RMSE()
                error_df['RMSE'] = rmse.evaluate(truth_mv_avg,pred_freq)



                
                score_df_list.append(error_df)
                
    score_df = pd.concat(score_df_list)



    #save score output to a csv file
    score_df.to_csv(f"../estimates/model_scores_output3.csv",index = False)


