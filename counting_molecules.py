#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 10:30:03 2018

@author: jack
"""

import pickle
import numpy as np

import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt

import nanonispy as nap

import mahotas

from scipy.spatial import distance
from scipy import optimize as _optimize

from sklearn.cluster import AffinityPropagation, Birch, AgglomerativeClustering

from skimage.draw import polygon
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import find_contours

import pairwise_chirality

### read sxm file, requires nanonispy

def read_data(filename, channel='Z', direction='forward'):
    if filename.endswith(".sxm"):
        scan = nap.read.Scan(filename)
        
        ## take the Z forward signal
        im = scan.signals[channel][direction]
        if scan.header['scan_dir'] == 'down':
            im = np.flipud(im)
    
        ## add mean values of image if scan is incomplete
        im[np.isnan(im)] = np.mean(im[~np.isnan(im)])
    
        ## factor for converting between pixels/ realspace
        rescale = scan.header['scan_range']/scan.header['scan_pixels']
    elif filename.endswith(".p"):
        d = pickle.load(open(filename, 'rb'))
        im = d['image']
        rescale = d['rescale']
        
    else:
        print('not an sxm file')
        return


    return im, rescale




### functions for doing a 2d plane fit to the image data:

def _plane(a0, a1, b1, x0, y0):
    return lambda x,y: a0 +a1*(x-x0) +b1*(y-y0)

def _planemoments(data):
    a0 = np.abs(data).min()
    index = (data-a0).argmin()
    x, y = data.shape
    x0 = float(index / x)
    y0 = float(index % y)
    a1 = 0.0
    b1 = 0.0
    return a0, a1, b1, x0, y0

def _fitplane(data):
    params = _planemoments(data)
    errorfunction = lambda p: np.ravel(_plane(*p)(*np.indices(data.shape)) - data)
    p, success = _optimize.leastsq(errorfunction, params)
    return p

def _return_plane(params, data):
    _fit_data = _plane(*params)
    return _fit_data(*np.indices(data.shape))

def _plane_fit_2d(scan_image):
    return scan_image - _return_plane(_fitplane(scan_image),scan_image)


def filter_image(im, gaussian_pixels=50):
    ## apply gaussian filtering
    im = im - gaussian(im, gaussian_pixels)
    ## normalize image
    im = im/np.amax(im)

    ## plane fit subtraction of image
    im = _plane_fit_2d(im)

    return im



### otsu to find the molecules in the filtered image

def get_contours(im, minimum_radius=.2e-9, minimum_separation=0, rescale=(1,1), zernike_radius=None):
    
    ## use otsu to find the threshold for molecules
    otsu_output = threshold_otsu(im)
    ## find the contours using the otsu threshold
    contours = find_contours(im, otsu_output)

    ##set contour length threshold:
    min_pixels = 2 * np.pi * np.sqrt(np.divide(minimum_radius**2, rescale[0]*rescale[1]))

    ## find contours larger than the minimum radius
    ## away from the edges of the image so they're "complete" molecules
    real_contours = []
    for c in contours:
        if min(c[:,0]) > 2 and max(c[:,0]) < im.shape[0]-2:
            if min(c[:,1]) > 2 and max(c[:,1]) < im.shape[1]-2:
                if len(c) > min_pixels:
                    real_contours.append(c)

    ## concatenate contours within minimimum_separation of each other
    new_contours = []
    used_indexes = []
    for ii, cc in enumerate(real_contours):
        for jj, cc2 in enumerate(real_contours):
            if jj>ii and ii not in used_indexes:
                dist = distance.cdist(cc*rescale, cc2*rescale)
                if np.amin(dist) < minimum_separation:
                    cc = np.concatenate((cc, cc2))
                    used_indexes.append(jj)
        if ii not in used_indexes:
            new_contours.append(cc)
            
    ### make a box that every contour will fit/ be rotatable in
    boxsize = 0
    diagonals = []
    for c in new_contours:
        x = c[:,1] #columns
        y = c[:,0] #rows
        xmin = int(np.min(x))
        xmax = int(np.max(x))
        ymin = int(np.min(y))
        ymax = int(np.max(y))
        if xmax - xmin > boxsize:
            boxsize = xmax - xmin
        if ymax - ymin > boxsize:
            boxsize = ymax - ymin
        diagonal = np.sqrt((xmin - xmax)**2 + (ymin - ymax)**2)
        diagonals.append(diagonal)
        if boxsize < diagonal:
            boxsize = int(diagonal)
    if boxsize/2 == 0:
        boxsize = boxsize + 1
    center = int(boxsize/2)
    
    
    #### create templates of all contour'd molecules; find contour length and maximum pixel height
    
    templates = []
    contour_lengths = []
    max_pixels = []
    
    for ii, c in enumerate(new_contours):
        poly = polygon(c[:,0], c[:,1])
        template = np.zeros((boxsize, boxsize))
        centerx_of_poly = int(np.mean(poly[0]))
        centery_of_poly = int(np.mean(poly[1]))
        translate_poly = (poly[0] - centerx_of_poly + center, poly[1] - centery_of_poly + center)
        template[translate_poly] = im[poly] - otsu_output
        templates.append(template)
        
        contour_lengths.append(len(c))
        max_pixels.append(np.amax(template) - otsu_output)
        
    
    ## normalize the contour lengths and the max pixel value
    contour_lengths = [xx / max(contour_lengths) for xx in contour_lengths]
    max_pixels = [xx / max(max_pixels) for xx in max_pixels]
    
    ### calculate the Zernike moments, add the normalized contour lengths and pixel_max to the moments arrays
    
    zernike_moments = []
    
    ### use median contour size as default zernike radius unless specified in function call
    if zernike_radius == None:
        zernike_radius = int(np.median(diagonals))
        

    for template, length, pixel in zip(templates, contour_lengths, max_pixels):
        answer = mahotas.features.zernike_moments(template, radius=zernike_radius)
        
        ## append the contour length and the maximum height into the zernike_moments
        
        answer = np.append(answer, pixel)
        answer = np.append(answer, length)
        
        zernike_moments.append(answer)
        
    ## convert list to np array
    zernike_moments = np.asarray(zernike_moments)
    

    return new_contours, otsu_output, templates, contour_lengths, max_pixels, zernike_moments

### use sklearn clustering to categorize the contours by clustering

def sort_contours(zernike_moments, damping=.7, exemplars=None, method=None, n_clusters=None):
    
    if exemplars == None and (method == None or method == 'Birch') and n_clusters == None:
        af = Birch(threshold=damping, n_clusters=n_clusters).fit(zernike_moments)
        return af.labels_
        
    elif n_clusters is not None and (method == None or method == 'AgglomerativeClustering'):
        af = AgglomerativeClustering(n_clusters=n_clusters).fit(zernike_moments)
        return af.labels_
        
    elif exemplars is not None or method == 'AffinityPropagation':
        preferences = -10 * np.ones(zernike_moments.shape[0])
        preferences[exemplars] = 0

        af = AffinityPropagation(damping=damping, preference=preferences).fit(zernike_moments)

        return af.labels_

def sort_chirality(templates, sorted_labels, nrotations=10, category_indexes=None):
    return pairwise_chirality.sort_chirality(templates, sorted_labels, nrotations=nrotations, category_indexes=category_indexes)

    
##### plotting functions

    
def make_fig(shape, dpi=96.):
    ''' return (fig, ax), without axes or white space '''
    h, w = shape
    dpi = float(dpi)
    fig = plt.figure()
    fig.set_size_inches(w/dpi, h/dpi, forward=True)
    ax = plt.Axes(fig, [0,0,1,1])
    ax.set_axis_off()
    fig.add_axes(ax)

    return fig, ax

#### plot grid of templates

def plot_template_grid(templates):
    normtemplates = []
    for temp in templates:
        vmax = np.max(temp)
        if vmax != 0:
            normtemplates.append(temp/vmax)
        else:
            normtemplates.append(temp)
    
    ntemplates = len(templates)
    griddim = int(np.ceil(np.sqrt(ntemplates)))
    empty = np.zeros(np.shape(templates[0]))
    for _ in range(griddim**2 - ntemplates):
        normtemplates.append(empty)
    
    
    ## plot an array of all the extracted templates
    
    pdim = griddim * np.shape(templates[0])[0]
    fig, ax = make_fig((pdim, pdim))
    chunks = [normtemplates[i:i+griddim] for i in range(0, len(templates), griddim)]
    templategrid = np.vstack([np.hstack(chunk) for chunk in chunks])
    plt.imshow(templategrid)
    
    return fig, ax

#### plot numbered contours

def plot_unsorted(im, real_contours, filename, rescale=(1,1)):
    plt.figure(figsize=(10,10))
    extent = (0, im.shape[0]*rescale[0], im.shape[1]*rescale[1], 0)
    plt.imshow(im, cmap='gray', extent=extent)
    plt.gca().set_xlabel('x (m)')
    plt.gca().set_ylabel('y (m)')
    for i, c in enumerate(real_contours):
        tempx = np.multiply(c[:,1], rescale[0])
        tempy = np.multiply(c[:,0], rescale[1])
        plt.plot(tempx, tempy, c='lime', linewidth=1)
        annotatex = np.mean(tempx)
        annotatey = np.mean(tempy)
        plt.text(annotatex, annotatey, i, color='red',
                    verticalalignment='center', horizontalalignment='center')
    temp = str(filename) + "_image_with_contours.png"
    print('Wrote file: {}'.format(temp))
    plt.savefig(temp, bbox_inches='tight')

    return


def plot_contours_histogram(im, contours, rescale, sorted_labels, saveplot=False, filename=None):

    partition = {k:0 for k in range(len(contours))}
    
    for index, ii in enumerate(sorted_labels):
        partition[index] = ii + 1
        
    ### make the plot
    
    plt.figure(figsize=(20,10))
    ax = plt.subplot(1,2,1)
    ax2 = plt.subplot(1,2,2)
    
    extent = (0, im.shape[0]*rescale[0], im.shape[1]*rescale[1], 0)
    ax.imshow(im, cmap='gray', extent=extent)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    
    cmap = matplotlib.cm.get_cmap('viridis')
    cmap2 = matplotlib.cm.get_cmap('plasma_r')
    

    numbins = max(partition.values())
    
    bins = np.bincount(list(partition.values()))
    newbins = sorted(bins, reverse=True)
    
    new_partition = {k:0 for k in range(len(contours))}
    next_identical = False
    for ii in range(len(newbins)):
        if next_identical == True:
            i = -1
            for j in range(2):
                i = list(bins).index(newbins[ii], i + 1)
            old_val = i
            next_identical = False
        else:
            old_val = list(bins).index(newbins[ii])
        if ii < len(newbins)-1:
            if newbins[ii+1] == newbins[ii]:
                next_identical = True
        mask = []
        for index, value in partition.items():
            if value == old_val:
                mask.append(index)
        for jj in mask:
            new_partition[jj] = ii
    partition = new_partition
    
    newbins = newbins[:numbins]
    
    colors = []
    for ii in range(numbins):
        if ii % 2 == 0:
            colors.append(cmap(ii/numbins))
        else:
            colors.append(cmap2(ii/numbins))
    
    for ii, c in zip(partition.keys(), contours):
        color = colors[partition[ii]]
        tempx = np.multiply(c[:,1], rescale[0])
        tempy = np.multiply(c[:,0], rescale[1])
        ax.plot(tempx, tempy, c=color, linewidth=1.5)
    
        #use this line to number all the molecules
        #ax.annotate(str(ii), xy=(tempx[0], tempy[0]), color='g')
    
    ax2.bar(range(len(newbins)), newbins, color=colors)
    
    # ax2.hist(bins, bins=len(bins))
    title = "total categories = " + str(max(partition.values()) + 1) + " total molecules = " + str(len(contours))
    ax2.set_title(title)
    ax2.set_xlabel('molecule category')
    ax2.set_ylabel('count')
    
    if saveplot == True:
        if filename is None:
            filename = 'output_histogram'
        savename = filename + '.png'
        plt.savefig(savename)
    return


#### interactive manual sorting

def manual_resorting(pickle_file=None, filename=None, im=None, contours=None, partition=None, manual_categories=None, rescale=(1,1)):
    if not pickle_file is None:
        d = pickle.load(open(pickle_file, 'rb'))
        filename = d['filename']
        im = d['image']
        contours = d['contours']
        partition = d['partition']
        rescale = d['rescale']
        
        if manual_categories is None:
            manual_categories = d['manual_categories']
    
    manual_categories = int(manual_categories)

    fig = plt.figure(figsize=(20,10))
    gs = matplotlib.gridspec.GridSpec(1,2)
    ax = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    extent = (0, im.shape[0]*rescale[0], im.shape[1]*rescale[1], 0)
    ax.imshow(im, cmap='gray', extent=extent) #, vmin=np.amin(im)*.1)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')

    cmap = matplotlib.cm.get_cmap('viridis')
    cmap2 = matplotlib.cm.get_cmap('plasma_r')


    if manual_categories is not None:
        numbins = int(manual_categories)
    else:
        numbins = max(partition.values())

    bins = np.bincount(list(partition.values()))
    newbins = sorted(bins, reverse=True)

    new_partition = {k:0 for k in range(len(contours))}
    next_identical = False
    for ii in range(len(newbins)):
        if next_identical == True:
            i = -1
            for j in range(2):
                i = list(bins).index(newbins[ii], i + 1)
            old_val = i
            next_identical = False
        else:
            old_val = list(bins).index(newbins[ii])
        if ii < len(newbins)-1:
            if newbins[ii+1] == newbins[ii]:
                next_identical = True
        mask = []
        for index, value in partition.items():
            if value == old_val:
                mask.append(index)
        for jj in mask:
            new_partition[jj] = ii
    partition = new_partition

    newbins = newbins[:numbins]
    errors = np.sqrt(newbins)

    colors = []
    for ii in range(numbins):
        if ii % 2 == 0:
            colors.append(cmap(ii/numbins))
        else:
            colors.append(cmap2(ii/numbins))
    colors.append('black')

    for ii, c in zip(partition.keys(), contours):
        if partition[ii] < len(colors):
            color = colors[partition[ii]]
        else:
            color = 'black'
        tempx = np.multiply(c[:,1], rescale[0])
        tempy = np.multiply(c[:,0], rescale[1])
        ax.plot(tempx, tempy, c=color, linewidth=1.5)

        #use this line to number all the molecules
        #ax.annotate(str(ii), xy=(tempx[0], tempy[0]), color='g')

        tempx = np.multiply(c[:,1], rescale[0])
        tempy = np.multiply(c[:,0], rescale[1])
        ax.plot(tempx, tempy, c=color, linewidth=1.5)

    bins = newbins
    bins = bins / np.sum(bins)
    errors = np.zeros(len(bins))

    rects = ax2.bar(range(len(bins)), bins, color=colors, yerr=errors)
    title = " total molecules = " + str(len(contours))
    ax2.set_title(title)
    ax2.set_xlabel('molecule category')
    ax2.set_ylabel('fraction of total count')
    
    for rect in rects:
        height = rect.get_height()
        height = height * len(contours)
        ax2.text(rect.get_x()+ rect.get_width()/2., 1.05*rect.get_height(), '%d' % int(height), ha='center', va='bottom')


    key_color = 0

    def onclick(event):
        global key_color
        xclick = event.xdata
        yclick = event.ydata
        ax.cla()
        ax.imshow(im, cmap='gray', extent=extent)#, vmin=np.amin(im)*.1)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')

        for ii, c in zip(partition.keys(), contours):
            xmin = min(c[:,1])*rescale[0]
            xmax = max(c[:,1])*rescale[0]
            ymin = min(c[:,0])*rescale[1]
            ymax = max(c[:,0])*rescale[1]

            if xclick < xmax and xclick > xmin and yclick < ymax and yclick > ymin:
                if not type(key_color) is int:
                    key_color = 0
                if key_color > len(colors):
                    color = 'black'
                else:
                    color = colors[key_color]
                partition[ii] = key_color
                tempx = np.multiply(c[:,1], rescale[0])
                tempy = np.multiply(c[:,0], rescale[1])
                ax.plot(tempx, tempy, c=color, linewidth=1.5)
            else:
                if partition[ii] < len(colors):
                    color = colors[partition[ii]]
                else:
                    color = 'black'
                tempx = np.multiply(c[:,1], rescale[0])
                tempy = np.multiply(c[:,0], rescale[1])
                ax.plot(tempx, tempy, c=color, linewidth=1.5)

        bins = np.bincount(list(partition.values()))
        errors = np.sqrt(bins)
        bins = bins / np.sum(bins)
        errors = np.zeros(len(bins))
        ax2.cla()
        rects = ax2.bar(range(len(bins)), bins, color=colors, yerr=errors)
        title = " total molecules = " + str(len(contours))
        ax2.set_title(title)
        ax2.set_xlabel('molecule category')
        ax2.set_ylabel('fraction of total count')
        
        for rect in rects:
            height = rect.get_height()
            height = height * len(contours)
            ax2.text(rect.get_x()+ rect.get_width()/2., 1.05*rect.get_height(), '%d' % int(height), ha='center', va='bottom')
        
        fig.canvas.draw()


    def keypress(event):
        global key_color, quit_fig
        try:
            key_color = int(event.key)
        except:
            key_color = 0
        if event.key == 'o':
            export = {'filename':filename, 'image':im, 'contours':contours, 'rescale':rescale, 'partition':partition, 'manual_categories':manual_categories}
            save_name = str(filename) + '_manually_sorted.p'
            f = open(save_name, 'wb')
            pickle.dump(export, f)
            f.close()
        if event.key == 'q':
            quit_fig = True


    fig.canvas.mpl_connect('key_press_event', keypress)
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    global quit_fig
    quit_fig = False    
    while quit_fig == False:
        plt.pause(.05)
        pass
    
    plt.close()
    
    return partition

### default sorting function

def default_sort(filename, sort_by_chirality=False):
    im, rescale = read_data(filename)
    im = filter_image(im)
    contours, otsu_output, templates, contour_lengths, max_pixels, zernike_moments = get_contours(im, rescale=rescale, minimum_separation=0)
    sorted_labels = sort_contours(zernike_moments, damping=.3, method='Birch')
    if sort_by_chirality == True:
        sorted_labels = sort_chirality(templates, sorted_labels)
    plot_contours_histogram(im, contours, rescale, sorted_labels, saveplot=True, filename=filename)
