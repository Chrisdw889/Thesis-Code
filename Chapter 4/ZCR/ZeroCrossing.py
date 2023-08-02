import os
import pickle
from multiprocessing import Pool

import numpy as np
import math

import pandas as pd

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

def event_search(traces, snr, thresh_val, window_size, step_size):
    '''
    ZCR detector main function.
    Locates events based on a precomputed threshold value.
    '''
    crossing_rates = calculate_ZCR(traces, window_size, step_size)
    evt_locs = []
    for rates in crossing_rates:
        evt_windows = one_window_search(rates, thresh_val, thresh_val)
        this_evt_locs = evt_windows_to_locs(evt_windows, window_size, step_size)
        evt_locs.append(this_evt_locs)
    return evt_locs

######## Background ########
def analyse_background(trace, window_size, step_size, bins=45):
    '''
    Fits a gaussian to the background ZCR distribution
    '''
    zcrs = calculate_ZCR(trace, window_size, step_size)
    bin_mids, counts, gauss_params, norm_factor = hist_analysis(zcrs, bins)
    
    return bin_mids, counts, gauss_params, norm_factor

def hist_analysis(zcrs, bins):
    '''
    analyse_background helper function
    '''
    flat_zcrs = np.array(flatten_zcrs_list(zcrs))
    
    counts, bin_edges = np.histogram(flat_zcrs, bins=bins)
    bin_mids = bin_edges[:-1] + (np.diff(bin_edges) / 2)
    
    params, cov = curve_fit(gauss, bin_mids, counts, p0=(8000, 0.5, flat_zcrs.std()))
    
    area = quad(gauss, -np.inf, np.inf, args=tuple(params))
    
    return bin_mids, counts, params, area[0]

def limit_finder(gauss_params, start, end, prob_threshold):
    '''
    Finds the lower integral limit that corresponds to a desired area given a precomputed distribution.
    '''
    trial_limits = np.linspace(start, end, 10)
    delta = abs(trial_limits[1] - trial_limits[0])
    
    areas = [quad(gauss, x, np.inf, args=tuple(gauss_params))[0] for x in trial_limits]
    
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

# One Window
def one_window_search(window_vals, threshold, return_threshold):
    '''
    Locates windows within an I(t) trace with anomalous values
    '''
    evt_locs = []
    
    i = 0
    while i < len(window_vals):
        val = window_vals[i]
        if val < threshold:
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
        
        if trace[i] > return_threshold:
            start = i + 1
            break
        elif i == 0:
            start = i
            break
        else:
            i -= 1
    
    j = end
    while j < max_idx:
        
        if trace[j] > return_threshold:
            end = j
            break
        elif j == (max_idx - 1):
            end = j + 1
            break
        else:
            j += 1
    
    return [start, end]
    
    
######## Utilities #########    

def gauss(x, A, mu, sig):
    y = A * np.exp((-1 * (x - mu) ** 2) / (2 * sig**2))
    return y

def flatten_zcrs_list(zcrs):
    '''
    Flattens a list of zcr arrays into a single axis.
    '''
    flat = []
    for entry in zcrs:
        flat.extend(entry)
    return np.array(flat)

def format_windows(data, window_size, step_size):
    '''
    Computes the possible windows for an I(t) trace given a window size and step size.
    '''
    windows = []
    
    i = 0
    j = i + window_size 
    
    while True: 
        if j > len(data):
            #windows.append(data[i:])
            break
        else:
            windows.append(data[i:j])
        
        i += step_size
        j = i + window_size
    return windows

def calculate_ZCR(data, window_size, step_size):
    '''
    Divides I(t) traces into possible windows and computes the ZCR for each.
    '''
    windows_all = []
    for trace in data:
        windows = format_windows(trace, window_size, step_size)
        windows_all.append(windows)
    windows_all = np.array(windows_all)
    
    zcrs = (np.diff(np.sign(windows_all), axis=2) != 0).sum(axis=2) / (window_size - 1)
    return zcrs
    

if __name__ == '__main__':
    pass