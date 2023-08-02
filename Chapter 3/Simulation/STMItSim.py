# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def STMItSim(event_freq, dur_mean, dur_sig, height_mean, height_sig, curr_sig_mean, curr_sig_sig, trace_num, points):
    '''
    Main function for constructing simulated I(t) data containing events. 
    '''
    
    # Noise event params
    noise_event_freq = 20
    noise_dur_mean = 2
    noise_dur_sig = 0.5
    noise_height_mean = 0.025
    noise_height_sig = 0.001
    noise_curr_sig_mean = 0.001
    noise_curr_sig_sig = 0
    
    
    traces = np.zeros((points, trace_num))
    trace_noise_locs = []
    trace_evt_locs = []
    for i in range(trace_num):
        arr = genBaseline(points)
        
        arr, noise_locs = genEvents(arr, noise_event_freq, noise_dur_mean, noise_dur_sig, noise_height_mean, noise_height_sig, noise_curr_sig_mean, noise_curr_sig_sig, points) # Noise
        arr, evt_locs = genEvents(arr, event_freq, dur_mean, dur_sig, height_mean, height_sig, curr_sig_mean, curr_sig_sig, points) # Molecule
        
        traces[:, i] = arr
        trace_noise_locs.append(noise_locs)
        trace_evt_locs.append(evt_locs)
        print(i)
        
    return traces, trace_noise_locs, trace_evt_locs
    

def genEvents(arr, event_freq, dur_mean, dur_sig, height_mean, height_sig, curr_sig_mean, curr_sig_sig, points):
    '''
    Populates a trace with I(t) events
    '''
    
    event_chance = event_freq / points
    rand = np.random.random(points)
    
    locations = []
    j = 0
    while j < points:
        if rand[j] < event_chance:
            event_current_std = np.random.normal(curr_sig_mean, curr_sig_sig)
            event_duration = int(np.random.normal(dur_mean, dur_sig))
            event_end = int(j + event_duration)
                
            if event_end > points:
                event_end = points
                event_duration = points - j
                    
            event_height = np.random.normal(height_mean, height_sig)
                
            event_currents = np.random.normal(event_height, event_current_std, event_duration)
                
            arr[j : event_end] = event_currents
            locations.append([j, event_end])
            j = event_end + 1
            continue
        j+=1
    return arr, locations
    
    
def genBaseline(points):
    '''
    Creates the baseline I(t) trace without molecular events or noise events.
    '''
     # Baseline point distribution
    setpoint = 0.0
    baseline_std = 0.005
    
    #Baseline sinusoid distribution
    inter_freq_mean = 500
    inter_freq_std = 100
    
    inter_amp_mean = 0.01
    inter_amp_std = 0.005
    
    #Frequency Modulation
    noise_freq_mean = 100
    noise_freq_std = 10 
    noise_freq_freq_mean = 10
    noise_freq_freq_std = 2
    
    noise_freq_fm_freq_mean = 10
    noise_freq_fm_freq_std = 0

    # Amplitude Modulation 
    noise_amp_freq_mean = 5
    noise_amp_freq_std = 1
    noise_amp_mean = 0.5
    noise_amp_std = 0.05
    
    noise_amp_fm_freq_mean = 10
    noise_amp_fm_freq_std = 2
    
    inter_amp = np.random.normal(inter_amp_mean, inter_amp_std)
    inter_freq = np.random.normal(inter_freq_mean, inter_freq_std)

    phase = np.random.uniform(0, 2 * np.pi)
    
    AM = genAmpMod(1, noise_amp_mean, noise_amp_std, noise_amp_freq_mean, noise_amp_freq_std, noise_amp_fm_freq_mean, noise_amp_fm_freq_std, points)
    FM = genAmpMod(1, noise_freq_mean, noise_freq_std, noise_freq_freq_mean, noise_freq_freq_std, noise_freq_fm_freq_mean, noise_freq_fm_freq_std, points)
  
    
    sin_one = [AM[i] * inter_amp * np.sin(FM[i] * inter_freq * x + phase) for i, x in enumerate(np.linspace(0, 1, points))]
    
    
    arr = np.random.normal(setpoint, baseline_std, points)
    
    arr = arr + sin_one
    
    plt.close('all')

    
    return arr

def genAmpMod(baseline, amp_mean, amp_std, FMAmp_mean, FMAmp_std, FMFreq_mean, FMFreq_std, points):
    '''
    Creates a continuous array of amplitude multipliers for modifying baseline traces.
    '''
    amp = np.random.normal(amp_mean, amp_std)
    FM = genFreqMod(FMAmp_mean, FMAmp_std, FMFreq_mean, FMFreq_std, points)
    phase = np.random.uniform(0, 2 * np.pi)
    
    AM = [baseline + (amp * np.sin(FM[i] * x + phase)) for i, x in enumerate(np.linspace(0, 1, points))]

    return AM

def genFreqMod(amp_mean, amp_std, FMFreq_mean, FMFreq_std, points):
    '''
    Creates a continuous array of frequency multipliers for modifying baseline traces.
    '''
    amp = np.random.normal(amp_mean, amp_std)
    FMFreq = np.random.normal(FMFreq_mean, FMFreq_std)
    phase = np.random.uniform(0, 2 * np.pi)
    
    FM = [amp * np.sin(FMFreq * x + phase) for x in np.linspace(0, 1, points)]
    
    return FM
    
if __name__ == '__main__':
    pass
    
