#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 11:31:15 2018

@author: jack
"""

import counting_molecules as ctm


filename = 'Ag111_APT_CO_044.sxm_output_dict.p'

## read in the image array and rescaling of pixels to distance
im, rescale = ctm.read_data(filename)

## Gaussian smooth and plane fit the image data:
im = ctm.filter_image(im)

## extract the contours, templates for each molecule; contour lenghts, maximum height and Zernike moments
contours, otsu_output, templates, contour_lengths, max_pixels, zernike_moments = ctm.get_contours(im, rescale=rescale, minimum_separation=0)

## plot the extracted contours, labelled by number, for selection of category exemplars
ctm.plot_unsorted(im, contours, filename, rescale=rescale)

## exemplar indices for the Helicene_Ag(111)008.sxm image
exemplars = [0, 1, 6, 10, 15, 36, 102]

## sort the molecules into categories
sorted_labels = ctm.sort_contours(zernike_moments, exemplars=exemplars)

#sorted_labels = ctm.sort_contours(zernike_moments, damping=.4, method='Birch', n_clusters=None)

## plot the sorted, histogram'd molecules:
ctm.plot_contours_histogram(im, contours, rescale, sorted_labels, saveplot=True, filename=filename)