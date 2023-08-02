import os
import pickle
from multiprocessing import Pool
from itertools import repeat

import numpy as np
import math

import pandas as pd

import tensorflow as tf
from tensorflow.keras.layers import Input, Conv1D, Dropout, Conv1DTranspose
from tensorflow.keras import Model
from tensorflow.keras.callbacks import EarlyStopping

from scipy.optimize import curve_fit
from scipy.integrate import quad

import matplotlib.pyplot as plt


def load_SNR_data(directory, files):
    '''
    Loads Simulation pkl files and extracts the raw I(t) traces and the event locations
    '''
    raw_traces = {}
    ground_truths = {} 
    
    for i, file in enumerate(files):
        print(i)
        with open(os.path.join(directory, file), 'rb') as f:
            data = pickle.load(f)
        SNR = float(file.split('_')[0])
        
        raw_traces[SNR] = data[0]
        ground_truths[SNR] = data[1]
        
    return raw_traces, ground_truths


def generate_model(window_size):
    '''
    Define the model architecture and optimizer
    '''
    inp =  Input((window_size, 1))

    x = Conv1D(32, 4, strides=1, activation='relu')(inp)
    x = Dropout(0.2)(x)
    # x = Conv1D(16, 7, padding='same', strides=2, activation='relu')(x)

    # y = Conv1DTranspose(16, 7, padding='same', strides=2, activation='relu')(x)
    # y = Dropout(0.1)(y)
    y = Conv1DTranspose(32, 4, strides=1, activation='relu')(x)

    out = Conv1DTranspose(1, 5, padding='same')(y)

    model = Model(inputs=inp, outputs=out)
    
    opt = tf.keras.optimizers.Adam(learning_rate=0.001)
    model.compile(optimizer=opt, loss='mse')
    
    return model

######## Background #########

def analyse_background(traces, model, window_size, step_size):
    '''
    Fits a gaussian to the background MAE distribution
    '''
    print('analysing')
    
    maes = calculate_maes(traces, model, window_size, step_size)
    bin_edges, counts, gauss_params, norm_factor = hist_analysis(maes)
    
    return maes, bin_edges, counts, gauss_params, norm_factor

def train_model(traces, window_size, step_size):
    '''
    Create and train a network model
    '''
    
    x_train = format_windows_all(traces, window_size, step_size)
    x_train = x_train.reshape((np.prod(x_train.shape[:-1]), x_train.shape[-1]))
    print(x_train.shape)
    
    model = generate_model(window_size)
    es = EarlyStopping(patience=1)
    
    history = model.fit(x_train, x_train, epochs=5, batch_size=128, validation_split=0.1, callbacks=[es])
    
    return model, history

def hist_analysis(maes):
    '''
    analyse_background helper function
    '''
    print(np.array(maes).shape)
    
    flat_maes = np.array(maes).flatten()
    
    counts, bin_edges = np.histogram(flat_maes, bins=100)
    bin_mids = bin_edges[:-1] + (np.diff(bin_edges) / 2)
    
    params, cov = curve_fit(gauss, bin_mids, counts, p0=(counts.max(), flat_maes.mean(), flat_maes.std()))
    
    area = quad(gauss, -np.inf, np.inf, args=tuple(params))
    
    return bin_edges, counts, params, area[0]


def limit_finder(gauss_params, start, end, prob_threshold):
    '''
    Finds the upper integral limit that corresponds to a desired area given a precomputed distribution.
    '''
    trial_limits = np.linspace(start, end, 10)
    delta = abs(trial_limits[1] - trial_limits[0])
    
    areas = [quad(gauss, -np.inf, x, args=tuple(gauss_params))[0] for x in trial_limits]
    
    dists = [abs(prob_threshold - x) for x in areas]
    
    closest_idx = np.argmin(dists)
    closest_area = areas[closest_idx]
    
    if  closest_area > (prob_threshold - 0.00005) and closest_area < (prob_threshold + 0.00005):
        return trial_limits[closest_idx]
    else:
        return limit_finder(gauss_params, 
                            trial_limits[closest_idx] - delta, 
                            trial_limits[closest_idx] + delta, 
                            prob_threshold)


######## Event Searching ########

def event_search(traces, model, thresh_val, window_size, step_size):
    '''
    AE detector main function.
    Locates events based on a precomputed threshold value.
    '''
    maes_all = calculate_maes(traces, model, window_size, step_size)
    
    evt_locs = []
    for maes in maes_all:
        maes_stepped = np.array(maes[::step_size])
        evt_windows = one_window_search(maes_stepped, thresh_val, thresh_val)
        this_evt_locs = evt_windows_to_locs(evt_windows, window_size, step_size)
        evt_locs.append(this_evt_locs)
    return evt_locs

def event_search_from_maes(maes_all, thresh_val, window_size, step_size, prestepped=False):
    '''
    AE detector alternate main function.
    Locates events based on a precomputed threshold value.
    Uses precomputed MAE values.
    '''
    
    evt_locs = []
    for maes in maes_all:
        if prestepped:
            maes_stepped = maes
        else:
            maes_stepped = np.array(maes[::step_size])
            
        evt_windows = one_window_search(maes_stepped, thresh_val, thresh_val)
        this_evt_locs = evt_windows_to_locs(evt_windows, window_size, step_size)
        evt_locs.append(this_evt_locs)
    return evt_locs


# One Window
def one_window_search(window_vals, threshold, return_threshold):
    '''
    Locates windows within an I(t) trace with anomalous values
    '''
    evt_locs = []
    
    i = 0
    while i < len(window_vals):
        val = window_vals[i]
        if val > threshold:
            locations = scope_event(window_vals, i, return_threshold)
            evt_locs.append(locations)
            i = locations[1]
        else:
            i += 1
            
    return np.array(evt_locs)

def scope_event(trace, idx, return_threshold):
    '''
    Locates the start and end of a flagged anomaly
    '''
    max_idx = len(trace)
    start, end = idx, idx + 1
    
    i = start
    while i >= 0:
        
        if trace[i] < return_threshold:
            start = i + 1
            break
        elif i == 0:
            start = i
            break
        else:
            i -= 1
    
    j = end
    while j < max_idx:
        
        if trace[j] < return_threshold:
            end = j
            break
        elif j == (max_idx - 1):
            end = j + 1
            break
        else:
            j += 1
    
    return [start, end]


def evt_windows_to_locs(evt_windows, window_size, step_size):
    '''
    Converts event locations from window index to point index
    '''
    evt_locs = []
    
    for start_window, end_window in evt_windows:
        start_idx = start_window * step_size
        end_idx = (end_window - 1) * step_size + window_size
        evt_locs.append([start_idx, end_idx])
        
    return evt_locs



######## Utilities #########

def gauss(x, A, mu, sig):
    y = A * np.exp((-1 * (x - mu) ** 2) / (2 * sig**2))
    return y


def format_windows(data, window_size, step_size):
    '''
    Computes the possible windows for an I(t) trace given a window size and step size.
    '''
    windows = []
    
    i = 0
    j = i + window_size 
    
    while True: 
        if j > len(data):
            break
        else:
            windows.append(data[i:j])
        
        i += step_size
        j = i + window_size
    return windows

def format_windows_all(traces, window_size, step_size):
    '''
    Computes the possible windows for a list of I(t) traces given a window size and step size.
    '''
    windows_all = []
    for trace in traces:
        windows = format_windows(trace, window_size, step_size)
        windows_all.append(windows)
    return np.array(windows_all)

def calculate_maes(traces, model, window_size, step_size):
    '''
    Divides I(t) traces into possible windows and computes the MAE for each.
    '''
    
    trace_maes = []
    for i, trace in enumerate(traces):
        trace_windows = format_windows(trace, window_size, step_size)
        
        maes = calculate_maes_helper(trace_windows, model)
        trace_maes.append(maes)
    
    return trace_maes

def calculate_maes_helper(windows, model):
    windows = np.array(windows)
    windows = np.expand_dims(windows, axis=2)
    recons = model.predict(np.array(windows))
    maes = np.mean(np.abs(recons - windows), axis=1)
    return maes

def locs_to_preds(evt_locs, traces):
    '''
    Converts event start and end locations to a one hot encoding
    '''
    predictions = np.zeros(traces.shape)
    
    for i, trace_evts in enumerate(evt_locs):
        for loc in trace_evts:
            start = loc[0]
            end = loc[1]
            
            predictions[i, start : end] = 1
            
    return predictions




