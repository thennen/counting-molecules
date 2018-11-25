#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 15:04:08 2018

@author: jack
"""

import counting_molecules as ctm


filename = 'Helicene_Ag(111)008.sxm'

## read in the image array and rescaling of pixels to distance
im, rescale = ctm.read_data(filename)

## Gaussian smooth and plane fit the image data:
im = ctm.filter_image(im)

## extract the contours, templates for each molecule; contour lenghts, maximum height and Zernike moments
contours_dict = ctm.get_contours(im, rescale=rescale)

## plot the extracted contours, labelled by number, for selection of category exemplars
ctm.plot_unsorted(im, contours_dict['contours'], filename, rescale=rescale)

## exemplar indices for the Helicene_Ag(111)008.sxm image
exemplars = [1, 2, 15, 21, 34, 52, 84]

## sort the molecules into categories
sorted_labels = ctm.sort_contours(contours_dict['zernike_moments'], exemplars=exemplars)

## plot the sorted, histogram'd molecules:
ctm.plot_contours_histogram(im, contours_dict['contours'], rescale, sorted_labels, saveplot=True, filename=filename)

## select the categories that need chiral sorting:
category_indexes=[0]

## re-sort the desired categories, separating by chirality
sorted_labels = ctm.sort_chirality(contours_dict['templates'], sorted_labels, category_indexes=category_indexes)

ctm.plot_contours_histogram(im, contours_dict['contours'], rescale, sorted_labels, saveplot=True, filename=filename)