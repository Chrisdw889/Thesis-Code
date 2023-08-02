import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

base_std = 1e-3 # Standard Deviation of I(t) traces baseline

def generate_traces(num_traces, length, SNR, event_dur, dur_dev, event_freq):
    '''
    Simulation main function.
    '''
    traces = []
    event_locs = []
    change_points = []
    
    for i in range(num_traces):
        trace, event_loc, change_point = generate_trace(length, SNR, event_dur, dur_dev, event_freq)
        traces.append(trace)
        event_locs.append(event_loc)
        change_points.append(change_point)
    
    return np.array(traces), np.array(event_locs), np.array(change_points)

def generate_trace(length, SNR, event_dur, dur_dev, event_freq):
    '''
    Simulation helper function.
    Generates a baseline then adds events.
    '''
    base = np.random.normal(0, base_std, length)
    base, locs, change_points = add_events(base, SNR, event_dur, dur_dev, event_freq)
    return base, locs, change_points

def add_events(base, SNR, event_dur, dur_dev, event_freq):
    '''
    Populate a baseline trace with events of user-defined SNR, duration, and frequency.
    '''
    evt_height = SNR * base_std
    
    event_starts = np.random.randint(0, len(base), event_freq)
    event_ends = []
    
    modifier_trace = np.zeros(len(base))
    
    for i in range(len(event_starts)):
        #generate end locations
        evt_dur = int(round(np.random.normal(event_dur, dur_dev, 1)[0]))
        event_end = event_starts[i] + evt_dur
        event_ends.append(event_end)
        
        #modify trace
        modifier_trace[event_starts[i] : event_end] = 1
        
    event_ends = np.array(event_ends)
    
    event_starts = np.expand_dims(event_starts, axis=1)
    event_ends = np.expand_dims(event_ends, axis=1)
    
    change_points = np.concatenate((event_starts, event_ends), axis=1)
    
    return base + (modifier_trace * evt_height), modifier_trace, change_points

if __name__ == '__main__':
    # Example of running the simulation and saving the data.
    SNR = 0.1
    
    data = generate_traces(10, 1000000, SNR, 10, 1, 20)
    
    with open('{}_SNR_data.pkl'.format(SNR), 'wb') as f:
        pkl.dump(data, f)
    
    
    