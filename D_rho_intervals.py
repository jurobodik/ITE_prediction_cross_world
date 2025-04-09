################################################################################
# D_rho_intervals.py
#
# Title: ITE Prediction Intervals under Cross-World Assumptions
# Author: jurajbodik.com
# Description:
#   This script provides an implementation of prediction intervals 
#   for Individual Treatment Effects (ITE) using the D_rho_intervals function 
#   under the Cross-World Assumption, characterized by the correlation parameter `rho`.
#
#   Core Function:
#     - `D_rho_intervals`:
#         * Estimates prediction intervals for ITEs at new covariate values
#         * Combines uncertainty from both treatment arms using a correlation-based formula
#         * Optionally adds bootstrap-based confidence intervals to account for CATE variability
#
#   Dependencies:
#     - numpy
#     - pandas
#     - lightgbm
#     - scikit-learn
#
#   Inputs:
#     - X: covariates (numpy array or pandas DataFrame), shape (n_samples, n_features)
#     - Y: observed outcomes (array-like), shape (n_samples,)
#     - treatment: treatment assignment (0/1), shape (n_samples,)
#     - new_points: covariates for which to predict ITEs (same format as X)
#     - rho: assumed correlation between potential outcomes (float in [-1, 1])
#
#   Optional Inputs:
#     - desired_coverage (float): target marginal coverage (default: 0.9)
#     - train_calib_split (float): proportion of data used for training vs calibration (default: 0.8)
#     - weighted_conformal (bool): whether to use propensity-weighted conformal inference (default: False)
#     - conformal (str): conformal method to use ('CQR' or 'naive') (default: 'CQR', see R version for PCS or CLEAR)
#     - CQR_qr (str): placeholder for QR engine; not used in current version (default: 'auto')
#     - CQR_type (str): conformity score type ('multiplicative', 'additive', or 'naive') (default: 'multiplicative')
#     - CATE_estimate (str): method for estimating the CATE ('T-learner') (default: 'T-learner', see R version for other options like causal_forest)
#     - add_confidence_intervals (bool): whether to compute bootstrap confidence intervals (default: True)
#     - number_of_bootstraps (int): number of bootstrap replicates (default: 50)
#
# Supported Features:
#     - Conformal prediction methods: CQR, Naive (PCS and CLEAR not yet implemented in Python)
#     - CATE estimation via T-learner (Causal Forest not included in this version)
#     - Bootstrap-based confidence intervals (optional but recommended)
#
#   Output:
#     - Dictionary with:
#         * 'new_points': the new covariate values (DataFrame)
#         * 'CATE': estimated conditional average treatment effects (pandas Series)
#         * 'lower': lower prediction bounds for ITEs (pandas Series)
#         * 'upper': upper prediction bounds for ITEs (pandas Series)
################################################################################

import numpy as np
import pandas as pd
import lightgbm as lgb #We use this for quantile regression: change 'fit_lgb_quantile' for different estimation
from sklearn.ensemble import GradientBoostingClassifier #We use this for propensity score estimation
import catboost as cb
from pygam import ExpectileGAM

def D_rho_intervals(
  X,                                 # Feature matrix (pd.DataFrame)
  Y,                                 # (factual) outcome variable (np.array)
  treatment,                         # Treatment assignment (0 or 1) (np.array)
  new_points,                        # New data points for prediction - where do you want to predict? (pd.DataFrame)
  rho,                               # Cross-world correlation parameter
  desired_coverage=0.9,              # Target marginal coverage level for the prediction intervals
  train_calib_split=0.8,             # Fraction of data used for training (rest used for calibration)
  weighted_conformal=False,          # Whether to use propensity-score weighted conformal prediction
  conformal='CQR',                   # Conformal method: 'CQR' or 'naive'
  CQR_qr='auto',                     # Quantile regression method: 'qgam', 'lightgbm', 'catboost', 'randomforest'. 'auto' chooses qgam for d<=5, lightgbm otherwise
  CQR_type='multiplicative',         # Type of conformity score: 'multiplicative' or 'additive'
  CATE_estimate='T-learner',         # Method to estimate CATE
  add_confidence_intervals=True,     # Whether to add bootstrap confidence intervals
  number_of_bootstraps=50,           # Number of bootstrap replicates for CI
):
    def fit_quantile(X_train, Y_train, alpha):
        X_train = np.asarray(X_train, dtype=np.float64)
        Y_train = np.asarray(Y_train).flatten()
        mask = np.isfinite(X_train).all(axis=1) & np.isfinite(Y_train)
        X_train, Y_train = X_train[mask], Y_train[mask]
        method = CQR_qr
        if CQR_qr == 'auto':
            method = 'qgam' if X_train.shape[1] <= 5 else 'lightgbm'
            
        if method == 'qgam':
            model = ExpectileGAM(expectile=alpha).gridsearch(X_train, Y_train)
            return lambda X: model.predict(np.asarray(X, dtype=np.float64))
        elif method == 'lightgbm':
            params = {'objective': 'quantile', 'alpha': alpha, 'verbosity': -1}
            dtrain = lgb.Dataset(X_train, label=Y_train)
            model = lgb.train(params, dtrain, num_boost_round=100)
            return lambda X: model.predict(np.asarray(X))
        elif method == 'catboost':
            model = cb.CatBoostRegressor(loss_function=f'Quantile:alpha={alpha}', iterations=100, verbose=False)
            model.fit(X_train, Y_train)
            return lambda X: model.predict(np.asarray(X))
        elif method == 'randomforest':
            model = RandomForestRegressor(n_estimators=100, min_samples_leaf=5)
            model.fit(X_train, Y_train)
            return lambda X: np.quantile([tree.predict(np.asarray(X)) for tree in model.estimators_], q=alpha, axis=0)
        else:
            raise ValueError("Unsupported quantile regression method")

    def rho_distance(x, y, rho):
        return np.sqrt(x**2 + y**2 - 2 * rho * x * y)

    def CI_addition(x, rho):
        return x * ((1 + rho) / 2) ** 2

    def predict_interval(X, Y, new_points):
        X, new_points = np.array(X), np.array(new_points)
        Y = np.array(Y).flatten()
        n = X.shape[0]
        train_size = int(n * train_calib_split)
        train_idx, calib_idx = np.arange(train_size), np.arange(train_size, n)
        X_train, Y_train = X[train_idx], Y[train_idx]
        X_calib, Y_calib = X[calib_idx], Y[calib_idx]
        all_X = np.vstack([X, new_points])

        pred_95 = fit_quantile(X_train, Y_train, (1 + desired_coverage) / 2)(all_X)
        pred_50 = fit_quantile(X_train, Y_train, 0.5)(all_X)
        pred_05 = fit_quantile(X_train, Y_train, (1 - desired_coverage) / 2)(all_X)

        pred_95 = np.maximum(pred_95, pred_50)
        pred_05 = np.minimum(pred_05, pred_50)

        if CQR_type == 'multiplicative':
            l, u = pred_50 - pred_05, pred_95 - pred_50
            scores = [max((pred_50[i] - Y_calib[j]) / l[i] if l[i] else 0,
                          (Y_calib[j] - pred_50[i]) / u[i] if u[i] else 0) for j, i in enumerate(calib_idx)]
        else:
            scores = np.abs(pred_50[calib_idx] - Y_calib)

        gamma = np.quantile(scores, desired_coverage)

        if CQR_type == 'multiplicative':
            lower = pred_50 - gamma * (pred_50 - pred_05)
            upper = pred_50 + gamma * (pred_95 - pred_50)
        else:
            lower = pred_05 - gamma
            upper = pred_95 + gamma

        m = new_points.shape[0]
        return {'hat_f': pred_50[-m:], 'lower': lower[-m:], 'upper': upper[-m:]}

    X, Y, treatment = np.array(X), np.array(Y), np.array(treatment)
    new_points = np.array(new_points)
    data = pd.DataFrame(X)
    data['Y'] = Y
    data['treatment'] = treatment
    data0, data1 = data[data.treatment == 0], data[data.treatment == 1]

    est0 = predict_interval(data0.iloc[:, :-2], data0['Y'], new_points)
    est1 = predict_interval(data1.iloc[:, :-2], data1['Y'], new_points)

    hat_f0, hat_f1 = est0['hat_f'], est1['hat_f']
    lower_u0, upper_u0 = hat_f0 - est0['lower'], est0['upper'] - hat_f0
    lower_u1, upper_u1 = hat_f1 - est1['lower'], est1['upper'] - hat_f1

    CATE = hat_f1 - hat_f0
    lower = CATE - rho_distance(lower_u0, upper_u1, rho)
    upper = CATE + rho_distance(upper_u0, lower_u1, rho)

    index_names = [f"new_point_{i+1}" for i in range(new_points.shape[0])]
    return {
        'new_points': pd.DataFrame(new_points, columns=[f"X{i+1}" for i in range(new_points.shape[1])]),
        'CATE': pd.Series(CATE, index=index_names),
        'lower': pd.Series(lower, index=index_names),
        'upper': pd.Series(upper, index=index_names)
    }
