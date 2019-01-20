#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 20:46:30 2018

@author: jack
"""

import numpy as _np

from skimage.transform import rotate as _imrotate

from scipy.interpolate import interp2d as _interp2d

from scipy.optimize import brute as _brute
from scipy.optimize import minimize as _minimize


def _template_stack(template, nrotations=10):
    """ Returns 3d array of nrotations images rotated through full circle. """
    angles = _np.linspace(0, 360, nrotations)
    
    ## imrotate does a weird re-scaling of the data, necessary for later comparison
    h,w = template.shape
    stack = _np.empty((nrotations, h, w))
    for ii, rr in enumerate(angles):
        stack[ii,:,:] = _imrotate(template, rr)
    
    return stack


def _template_compare(template1, template2, x, y, r):
    """ Returns normalized pixel-to-pixel comparison of two image arrays with relative angle r. """
    _, h, w = template2.shape
    I2d = _interp2d(range(h), range(w), template2[0,:,:])
    t1 = _imrotate(template1[0,:,:], r)
    t2 = I2d(_np.arange(h) - y, _np.arange(w) - x)
    return 1 - _np.sum(_np.abs(t1 - t2))/_np.sum(_np.abs(t1)+_np.abs(t2))

def _template_compare_fixed_r(template1, template2, x, y, r_index):
    """  Returns normalized pixel-to-pixel comparison of two image arrays with fixed relative angle r_index. """
    nrotations = template1.shape[0]
    r_index = int(r_index)
    if r_index >= nrotations:
        r_index = nrotations - 1
    _, h, w = template2.shape
    I2d = _interp2d(range(h), range(w), template2[0,:,:])
    t1 = template1[r_index,:,:]
    t2 = I2d(_np.arange(h) - y, _np.arange(w) - x)
    return 1 - _np.sum(_np.abs(t1 - t2))/_np.sum(_np.abs(t1)+_np.abs(t2))


def _compare_templates(template, template2):
    """ Returns highest correlation between two pre-rotated image stacks. """
    nrotations = template.shape[0]

    ## optimization with imrotate
    def opt_func(input):
        x = input[0]
        y = input[1]
        r = input[2]
        return 1 - _template_compare(template, template2, x, y, r)
    
    ## optimization with pre-computed rotation stack
    def brute_opt(input):
        x = input[0]
        y = input[1]
        r = input[2]
        return 1 - _template_compare_fixed_r(template, template2, x, y, r)
    
    ranges = (slice(-4, 4, 2), slice(-4, 4, 2), slice(0, nrotations-1, 1))
    answer = _brute(brute_opt, ranges=ranges)
    answer[2] = answer[2]/nrotations*360
    answer = _minimize(opt_func, x0=answer)
    answer = answer['x']
    corr = 1 - opt_func(answer)
    
    return corr


def sort_chirality(templates, sorted_labels, nrotations=10, category_indexes=None):
    """ Returns list of chiral sorted contour index labels. """
    new_labels = [ii for ii in sorted_labels]
    templates = _np.array(templates)
    
    if category_indexes == None:
        category_indexes = [ii for ii in range(max(sorted_labels))]
    
    for category in category_indexes:
    
        mask = _np.array(sorted_labels) == category
        subset = templates[mask == True]
        
        ## pre-rotate all the templates
        subset = [_template_stack(ii, nrotations=nrotations) for ii in subset]
        
        correlations = _np.zeros(len(subset))
    
        for ii, template in enumerate(subset):
            original = _compare_templates(subset[0], template)
            mirror = _compare_templates(subset[0], _np.fliplr(template))
    
            if mirror > original:
                correlations[ii] = 1
    
        new_index = category + max(sorted_labels) + 1   
        jj = 0
        for ii, mask_val in enumerate(mask):
            if mask_val == True:
                if correlations[jj] == 1:
                    new_labels[ii] = new_index
                jj = jj + 1
    
    return new_labels
