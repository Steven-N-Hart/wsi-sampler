import numpy as np
import openslide
import skimage.morphology as morph
from scipy import ndimage
from skimage import color
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops


def getScaledCoordinates(xy, hw_in, hw_out):
    '''
    Compute the coordinates corresponding to a scaled version of an image

    Parameters
    ==========
    xy:
        list of (x,y) coordinates
    hw_in:
        [h_in, w_in] shape of the image to which the input coordinates refer
    hw_out:
        [h_out, w_out] shape of the scaled image where to compute the output coordinates

    Return
    ======
    xy_out:
        list of (x,y) coordinates in the output image
    '''
    assert len(hw_in) == 2
    assert len(hw_out) == 2
    assert len(xy) >= 1

    h_in = hw_in[0]
    w_in = hw_in[1]
    h_out = hw_out[0]
    w_out = hw_out[1]

    xy_out = []
    for XY in xy:
        x_out = int(np.floor(h_out * XY[0] / h_in))
        y_out = int(np.floor(w_out * XY[1] / w_in))
        xy_out.append([x_out, y_out])
    return (xy_out)


# select tissue region
def getTissueRegion(image, thumbnail_size=1000, return_thumbnail=False):
    '''
    Compute the box containing the tissue

    Parameters
    ==========
    image:
        OpenSlide image
    thumbnail_size:
        int, size of the thumbnail where to compute the box
    return_thumbnail:
        boolean, whether to return the thumbnail image


    Return
    ======
    box_coords:
        [x_ul, y_ul, x_br, y_br] coordinates of the box containing the tissue
    thumbnail:
        PIL.Image of the tissue in the thumbnail (if return_thumbnail=True)
    '''
    assert isinstance(image, openslide.OpenSlide), "input image should be an openslide wsi"

    w_out, h_out = image.dimensions  # ! openslide image dimensions: WxH

    thumb = np.array(image.get_thumbnail((thumbnail_size, thumbnail_size)))
    h_in, w_in, ch = thumb.shape

    thumb = color.rgb2gray(thumb)
    TH = threshold_otsu(thumb)
    thumb_BW = thumb < TH
    Strel = morph.disk(3)
    thumb_BW_dilated = morph.dilation(thumb_BW, Strel)
    thumb_BW_filled = ndimage.binary_fill_holes(thumb_BW_dilated, structure=np.ones((5, 5))).astype(int)

    # get biggest region
    L = label(thumb_BW_filled)
    props = regionprops(L)

    areas = []
    boxes = []
    centers = []
    for prop in props:
        areas.append(prop.area)
        boxes.append(prop.bbox)
        centers.append(prop.centroid)

    biggest_region = np.argmax(areas)
    x_ul, y_ul, x_br, y_br = boxes[biggest_region]

    out_coords = getScaledCoordinates([[x_ul, y_ul], [x_br, y_br]], [h_in, w_in], [h_out, w_out])
    out_coords = [out_coords[0][0], out_coords[0][1], out_coords[1][0], out_coords[1][1]]
    if return_thumbnail:
        thumb_out = thumb[x_ul:x_br, y_ul:y_br]
        return (out_coords, thumb_out)
    else:
        return (out_coords)


if __name__ == '__main__':
    SVS = '/mnt/d/TCGA/78c0e50f-116d-4b2d-bad0-89385ada6e6f/TCGA-A8-A09Z-01Z-00-DX1.D56497BE-5099-4537-9600-60F3213F7BF5.svs'
    img = openslide.OpenSlide(SVS)  # open the image
    box_tissue, thumbnail = getTissueRegion(img, return_thumbnail=True)  # get the tissue box, scaled level0
    thumbnail.show()
