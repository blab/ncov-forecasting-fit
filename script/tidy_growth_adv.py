from statistics import mode
import pandas as pd
import os
import numpy as np
import itertools




if __name__=='__main__':
    locations = ["USA","Japan", "United Kingdom"]
    models = ["GARW", "MLR", "FGA", "Piantham"]
    dates = ['2022-04-15','2022-04-22','2022-04-29','2022-05-06',
         '2022-05-13','2022-05-20','2022-05-27','2022-06-03',
         '2022-06-10','2022-06-17','2022-06-24','2022-06-30']



    #full model output set dict
    final_sets = {}

    #loop thorough different files of model versions
    for model in models:

        for location in locations:
            ga_dic = {}
            #filtering final_truth dataset to run location
            for date in dates:
                
                filepath = f"../plot-est/cast_estimates_full_{model}/{location}/full_growth_advantages_{date}.csv"
                
                #Check if file exists and continue if not
                if not os.path.exists(filepath):
                    continue

                #read models and add to dict
                ga_data = pd.read_csv(filepath)

                ga_dic[date] =  ga_data


            #dict of all ga datasets by pivot date
            final_sets[location, model] = ga_dic

            #print(final_sets)

            



    ga_df_list = []

    for model in models:

        for location in locations:

            for k, v in final_sets[location, model].items():
                pivot_date = pd.to_datetime(k)
                variants = final_sets[location, model][k]['variant']

                #creating ga dict
                ga_df = pd.DataFrame({'location': location, 'model': model, 'pivot_date': pivot_date, 'variant':variants})

                #unpacking prepped_data values
                if model =='GARW':
                    ga_df['date'] = v['date']
                ga_df[['median_ga', "ga_upper_80", "ga_lower_80" ]] = v[['median_ga', "ga_upper_80", "ga_lower_80" ]]

                
                ga_df_list.append(ga_df)




    #CONCAT ALL DATASETS
    ga_df = pd.concat(ga_df_list)


    #save score output to a csv file
    #output csv of all datasets for each country
    ga_df.to_csv(f"../estimates/ga_output.csv",index = False)


