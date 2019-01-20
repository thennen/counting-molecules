#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 11:31:15 2018

@author: jack
"""

import counting_molecules as ctm


filename = '111_test_data.p'

## read in the image array and rescaling of pixels to distance
im, rescale = ctm.read_data(filename)

## Gaussian smooth and plane fit the image data:
im = ctm.filter_image(im)

## extract the contours, templates for each molecule; contour lenghts, maximum height and Zernike moments
contours_dict = ctm.get_contours(im, rescale=rescale, minimum_separation=.5e-9, block_size=35, offset=-.5)

## plot the extracted contours, labelled by number, for selection of category exemplars
ctm.plot_unsorted(im, contours_dict['contours'], filename, rescale=rescale)

## exemplar indices for the 111 test data APT image
exemplars = [17, 39, 21, 5, 34, 9, 11, 47, 57]

## sort the molecules into categories
sorted_labels = ctm.sort_contours(contours_dict['zernike_moments'], exemplars=exemplars)

## less good sorting alternatives, no exemplar pre-selection required:
#sorted_labels = ctm.sort_contours(contours_dict['zernike_moments'], damping=.3, n_clusters=8)
#sorted_labels = ctm.sort_contours(contours_dict['zernike_moments'], damping=.35, method='Birch')

### plot the sorted, histogram'd molecules:
ctm.plot_contours_histogram(im, contours_dict['contours'], rescale, sorted_labels, saveplot=True, filename=filename)