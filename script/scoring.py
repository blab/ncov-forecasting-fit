from typing import final
import pandas as pd
import sklearn
from sklearn.metrics import log_loss, brier_score_loss
import os
import numpy as np
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod




def prep_freq_data(final_set):
    
    
    #return truth values as np.arrays
    
    truth_values = final_set['truth_freq'].to_numpy()
    
    #return nowcast and forecast values as np.arrays
    
    nowcast_values = final_set['median_freq_nowcast'].to_numpy()
    
    forecast_values = final_set['median_freq_forecast'].to_numpy()
    
    return truth_values, nowcast_values, forecast_values



class Scores(ABC):
    @abstractmethod
    def __init__(self):
        pass
class MAE(Scores):
    def __init__(self):
        pass
    def evaluate(self,truth, prediction):
        abs_error = np.abs(truth - prediction)
        return np.nanmean(abs_error)
        
class MSE(Scores):
    def __init__(self):
        pass
    def evaluate(self, truth, prediction):
        squared_error = np.square(truth - prediction)
        return np.nanmean(squared_error)
        
class LogLoss(Scores): #to-do
    def __init__(self):
        pass
        #mlr log loss error
    def evaluate(self, truth, prediction):
        log_loss = sklearn.metrics.log_loss(self.truth, prediction, eps=1e-15, normalize=True, sample_weight=None, labels=None)
        return log_loss


if __name__=='__main__':
    locations = ["USA","Japan"]
    models = ["GARW", "FGA"]
    dates = ['2022-01-24', '2022-02-04','2022-02-08','2022-02-18','2022-02-23',
         '2022-02-28','2022-03-03','2022-03-08','2022-03-15',
         '2022-03-21','2022-03-25','2022-04-07','2022-04-14','2022-04-27'
         ,'2022-05-06','2022-05-17','2022-05-20','2022-05-28','2022-06-09'
         ,'2022-06-14','2022-06-22']
    #Latest model run "truth"
    truth_set = pd.read_csv("omicron-countries-split_freq-combined-GARW.tsv", sep="\t")
    final_truth = truth_set.rename(columns = {'median_freq':'truth_freq'}, inplace = False)

    
    #full model output set dict


    final_sets = {}

    #loop thorough different files of model versions
    for model in models:

        for location in locations:
            pred_dic = {}
            #filtering final_truth dataset to run location
            location_truth = final_truth[final_truth['location']==location]
            location_truth = location_truth[['date','location','variant','truth_freq']]

            for date in dates:
                
                filepath = f"../estimates/cast_estimates_full_{model}/{location}/freq_full_{date}.csv"
                    
                #Check if file exists and continue if not
                if not os.path.exists(filepath):
                    continue
                    #read models and add to dict
                pred_dic[date] = pd.read_csv(filepath)
                
                    #loop through data and merge to final set
            #print(pred_dic['2022-02-28']['median_freq_nowcast'])
            #print(model,location)
            #print(location_truth['date'].min())
            #print(pred_dic['2022-02-28']['date'])
            final_sets_location = {k: pd.merge(location_truth,d) for k,d in pred_dic.items()}
            #print(final_sets_location['2022-02-28']['median_freq_nowcast'])
            final_sets[location, model] = final_sets_location
        
    #print(final_sets.keys())  
    #print(final_sets['USA','GARW']['2022-02-28']['median_freq_nowcast'])


    error_id_dict = {}
    #prep arrays of data

    for model in models:

        for location in locations:
            prepped_data = {k: prep_freq_data(v) for k,v in final_sets[location, model].items()}
              
            error_id_location = {}    
            for k, v in prepped_data.items():
                error_dict={}
                #MAE
                
                error = MAE()  
                error_dict['MAE'] = error.evaluate(v[0],v[1])
                #print(v[1])
                #MSE
                mse = MSE()
                error_dict['MSE'] = mse.evaluate(v[0],v[1])
                
                error_id_location[k] = error_dict
            error_id_dict[location, model] = error_id_location
        
    #print(prepped_data)
    #print(error_id_dict)
    
#formatting output to desired format
    
#convert output to a dataframe
score_df = pd.DataFrame.from_dict({(location,model): error_id_dict[location][model] 
for location in error_id_dict.keys() 
for model in error_id_dict[location]}, orient='index')


def clean_score(score_df):
    #transform dict to columns
    scores = score_df.explode('MSE').reset_index(drop=False)
    scores.columns = ['location','date','Mean Abs Error','Mean Sq Error']
    #unnesting location and model columns
    iden_model = pd.DataFrame(scores["location"].to_list(), columns=['locations', 'model'])
    scores_final = pd.concat([iden_model,scores], axis = 1).drop(['location'], axis = 1)

    #save score output to a csv file
    scores_final.to_csv(f"../estimates/model_scores_outputt.csv",index = False)

clean_score(score_df)