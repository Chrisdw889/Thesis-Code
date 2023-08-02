import os
import pickle

import numpy as np

import pandas as pd


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

def event_search(background_trace, traces, st_dev, snr):
    '''
    Amplitude thresholding main function
    Calculates background average and standard deviation for deciding threshold
    Then locates events
    '''
    std = np.std(background_trace.flatten())
    avg = np.mean(background_trace.flatten())
    
    threshold = avg + (st_dev * std)
    
    evt_all_locations = []
    # For each trace
    for idx, trace in enumerate(traces):
        print('Trace index: ' + str(idx))
        evt_locations = find_events(trace, threshold, avg, st_dev, snr, idx)
        evt_all_locations.append(evt_locations)
    
    return evt_all_locations

def find_events(trace, threshold, return_threshold, st_dev, snr, idx):
    '''
    Locates events within an I(t) trace
    '''
    
    evt_locs = []
    
    i = 0
    while i < len(trace):
        val = trace[i]
        #if i > 900000 and i != 999999:
        #print(i)
        if val > threshold:
            locations = scope_event(trace, i, return_threshold)
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
            
            
if __name__ == '__main__':
    pass
        
    
    
        
    
    
    
    
    
        
        
            
        