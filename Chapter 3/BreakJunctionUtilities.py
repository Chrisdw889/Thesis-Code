# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:23:37 2022

@author: cxw416-admin
"""
import scipy.io as spio
import numpy as np

def formatData(t,s):
    if not isinstance(t,dict) and not isinstance(t,list):
        print("\t"*s+str(t))
    else:
        for key in t:
            print("\t"*s+str(key))
            if not isinstance(t,list):
                formatData(t[key],s+1)

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it fixes the issue of not properly recovering python dictionaries
    from mat files. It calls the function check keys to fix all entries
    which are still mat-objects.
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries.
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
        elif isinstance(dict[key], np.ndarray):
            dict[key] = _check_list(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries.
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        elif isinstance(elem, np.ndarray):
            dict[strg] = _check_list(elem)
        else:
            dict[strg] = elem
    return dict

def _check_list(list_obj):
    '''
    converts lists of mat-objects to lists of dictionaries.
    '''
    new_list = []
    for el in list_obj:
        if isinstance(el, spio.matlab.mio5_params.mat_struct):
            new_list.append(_todict(el))
        else:
            new_list.append(el)
    return new_list