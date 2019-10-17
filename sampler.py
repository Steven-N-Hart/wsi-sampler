import glob
import math
import os

import numpy as np
import openslide
from skimage.filters import threshold_otsu

#
# -----------------------------------------------------------
# Create function
def get_thumbnail(img, x_level0, y_level0, patch_size, output_path):
    """

    :param img:
    :param x_level0:
    :param y_level0:
    :param patch_size:
    :param output_path:
    :return:
    """
    patch = img.read_region((x_level0, y_level0), 0, (patch_size, patch_size))
    patch = patch.convert('RGB')
    fname = img.properties['aperio.Filename'].replace(' ', '_')
    fname += '_' + str(x_level0)
    fname += '_' + str(y_level0)
    fname += '_' + '0'
    fname += '_' + str(patch_size)
    fname += '.jpg'
    patch_arr = np.array(patch.convert('L'))
    # print('Getting close: {}'.format(np.average(patch_arr)))
    if np.average(patch_arr) < 230:
        global num_patches
        num_patches += 1
        patch.save(os.path.join(output_path, fname))


def process_svs(SVS,
                normalization_factor=1000,
                patch_size=512,
                buffer=10,
                output_path='/projects/shart/digital_pathology/data/TCGA/test_svs/test/wsi_sampler'):
    img = openslide.OpenSlide(SVS)

    global num_patches
    num_patches = 0

    # Normalize the thumbnail to the actual slide size
    x_max, y_max = img.dimensions
    thumbnail = img.get_thumbnail((normalization_factor, normalization_factor))

    grey_thumbnail = np.array(thumbnail.convert("L"))
    thresh = threshold_otsu(grey_thumbnail)
    mask = np.array(grey_thumbnail) < thresh

    # how many pixels in the raw image per pixel in mask
    x_num_orgPix_per_thumbPix = math.ceil(x_max / mask.shape[0])
    y_num_orgPix_per_thumbPix = math.ceil(y_max / mask.shape[1])
    # print(x_num_orgPix_per_thumbPix, y_num_orgPix_per_thumbPix)

    # Find out how many pixels in image mask to count as a patch in original
    num_x_mask_pixels_per_rawPatch = math.ceil(patch_size / x_num_orgPix_per_thumbPix)
    num_y_mask_pixels_per_rawPatch = math.ceil(patch_size / y_num_orgPix_per_thumbPix)
    # print(num_x_mask_pixels_per_rawPatch, num_y_mask_pixels_per_rawPatch)
    #Naresh: Adjust the buffer size to the Mask size ratio
    buffer_x = math.ceil(buffer / x_num_orgPix_per_thumbPix)
    buffer_y = math.ceil(buffer / y_num_orgPix_per_thumbPix)
    mask_x, mask_y = mask.shape
    x_mask_prev = 0

    # Iterate through the mask to identify positive pixels
    #Naresh:Adjusted the step size  to Mask size
    for x in range(buffer_x, mask_x - (buffer_x), num_x_mask_pixels_per_rawPatch):
        x_mask_window = x + num_x_mask_pixels_per_rawPatch
        if x_mask_window <= x_mask_prev:
            continue
        y_mask_prev = 0
        # Naresh:Adjusted the step size  to Mask size
        for y in range(buffer_y, mask_y - (buffer_y), num_y_mask_pixels_per_rawPatch):
            y_mask_window = y + num_y_mask_pixels_per_rawPatch
            # print('Evaluate: {} {} & {}'.format(y, y_mask_window, y_mask_prev))
            if y_mask_window <= y_mask_prev:
                continue
            if y % 100 == 0:
                print('X: {}\tY:{} of {} with total of {} so far'.format(x, y, mask.shape, num_patches), end='\r',
                      flush=True)
            if np.sum(mask[x:x_mask_window, y:y_mask_window]) > 0:
                # convert mask coordinates to level0 coordinates
                x_level0 = x * x_num_orgPix_per_thumbPix
                y_level0 = y * y_num_orgPix_per_thumbPix
                get_thumbnail(img, x_level0, y_level0, patch_size, output_path)
            # print('yamsk windoe: {}'.format(y_mask_window))
            y_mask_prev = y_mask_window
        x_mask_prev = x_mask_window

    print('Printed {} from {}'.format(num_patches, SVS))


num_patches = 0

if __name__ == '__main__':
    for SVS in glob.glob('/projects/shart/digital_pathology/data/TCGA/test_svs/test/*svs'):
        process_svs(SVS)
