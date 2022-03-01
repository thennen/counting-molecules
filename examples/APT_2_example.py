#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 11:31:15 2018

@author: jack
"""

import counting_molecules as ctm


filename = 'Ag111_APT_CO_044.p'

## read in the image array and rescaling of pixels to distance
im, rescale = ctm.read_data(filename)

## Gaussian smooth and plane fit the image data:
im = ctm.filter_image(im)

## extract the contours, templates for each molecule; contour lengths, maximum height and Zernike moments
contours_dict = ctm.get_contours(im, rescale=rescale, minimum_separation=.5e-9, block_size=35, offset=-.5)

## plot the extracted contours, labelled by number, for selection of category exemplars
ctm.plot_unsorted(im, contours_dict['contours'], filename, rescale=rescale)

## exemplar indices for the APT 044 image
exemplars = [28, 29, 5, 16, 13, 2]

## sort the molecules into categories
sorted_labels = ctm.sort_contours(contours_dict['zernike_moments'], exemplars=exemplars)

## plot the sorted, histogram'd molecules:
ctm.plot_contours_histogram(im, contours_dict['contours'], rescale, sorted_labels, saveplot=True, filename=filename)