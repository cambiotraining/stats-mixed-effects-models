import os
import glob
from datetime import datetime
from typing import Type
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from plotnine import *
import patchworklib as pw

def dgplotsknitr(results: Type[sm.regression.linear_model.RegressionResultsWrapper]) -> None:
    if isinstance(results, sm.regression.linear_model.RegressionResultsWrapper) is False:
        raise TypeError("Please provide a model fit.")
    else:
        
        residuals = results.resid.rename("residuals")
        predicted_values = results.fittedvalues.rename("predicted_values")
        std_resid = pd.Series(np.sqrt(np.abs(results.get_influence().resid_studentized_internal))).rename("std_resid")
        influence = results.get_influence()
        cooks_d = pd.Series(influence.cooks_distance[0]).rename("cooks_d")
        leverage = pd.Series(influence.hat_matrix_diag).rename("leverage")
        obs = pd.Series(range(len(residuals))).rename("obs")
        n_obs = len(obs.index)
        
        # combine Series into DataFrame
        model_values = residuals.to_frame().join(predicted_values).join(std_resid).join(cooks_d).join(leverage).join(obs)
        # add the total number of observations
        model_values["n_obs"] = n_obs
        
        p1 = (
        ggplot(model_values, aes(x = "predicted_values", y = "residuals"))
        + geom_point()
        + geom_smooth(se = False, colour = "red")
        + labs(title = "Residuals plot")
        + xlab("predicted values")
        + ylab("residuals")
        + theme_bw()
        )

        p2 = (
        ggplot(model_values, aes(sample = "residuals"))
        + stat_qq()
        + stat_qq_line(colour = "blue")
        + labs(title = "Q-Q plot")
        + xlab("theoretical quantiles")
        + ylab("sample quantiles")
        + theme_bw()
        )

        p3 = (
        ggplot(model_values, aes(x = "predicted_values", y = "std_resid"))
        + geom_point()
        + geom_smooth(se = False, colour = "red")
        + labs(title = "Location-Scale plot")
        + xlab("predicted values")
        + ylab(u"\u221A"'|standardised residuals|')
        + theme_bw()
        )

        p4 = (
        ggplot(model_values, aes(x = "obs", y = "cooks_d"))
        + geom_point()
        + geom_segment(aes(xend = "obs", yend = 0), colour = "blue")
        + geom_hline(aes(yintercept = 0))
        + geom_hline(aes(yintercept = 4/n_obs), colour = "blue", linetype = "dashed")
        + labs(title = "Influential points")
        + xlab("observation")
        + ylab('cook\'s d')
        + theme_bw()
        )

        p1 = pw.load_ggplot(p1, figsize=(3,2))
        p2 = pw.load_ggplot(p2, figsize=(3,2))
        p3 = pw.load_ggplot(p3, figsize=(3,2))
        p4 = pw.load_ggplot(p4, figsize=(3,2))

        dgplots = (p1 | p2) / (p3 | p4)
        date = datetime.now().strftime("%Y_%m_%d-%I-%M-%S_%p")
        dgplots.savefig(fname=(f"images/dgplots/{date}_dgplots" ".png"))
        # get list of images created by dgplots
        list_of_files = glob.glob('images/dgplots/*dgplots.png')
        # get the latest file
        latest_file = max(list_of_files, key=os.path.getctime)
        
        return latest_file
