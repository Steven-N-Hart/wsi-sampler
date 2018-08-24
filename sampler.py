import openslide
import numpy as np


class WsiSample(object):

    """This class will scan the highest level of a whole slide image to identify regions with sufficient tissue
    Args:
        wsi_path: '/path/to/wsi.svs' [Required]
        threshold = Value at which one considers empty space [200]
        x_buffer = How many pixels to ignore along the x-axis [50]
        y_buffer = How many pixels to ignore along the y-axis [50]
        patch_buffer = How many pixels to include before the tissue is detected [5]
        patch_size = How large of a coordinate set to return [300]
        patch_interval = Fraction for how much the patches should overlap [0.1]


    Usage:
        import WsiSample
        import openslide

        wsi_path='/path/2/wsi.svs'
        patch_size = 300
        level = 0

        obj = WsiSample(wsi_path='/path/2/wsi.svs', patch_size=patch_size)
        OSobj = openslide.OpenSlide(wsi_path)

        for x,y in obj.get_targets(): # x and y are coordinates to extract from your OpenSlide Object
            img_patch = OSobj.read_region((x,y), level, (patch_size, patch_size))
    """

    def __init__(self,
                 wsi_path=None,
                 patch_interval=0.1,
                 threshold=200,
                 x_buffer=50,
                 y_buffer=50,
                 patch_size=300,
                 patch_buffer=5,
                 level=None
                 ):
        self.wsi_path = wsi_path  # '/path/to/wsi.svs'
        self.threshold = threshold  # Value at which one considers empty space
        self.x_buffer = x_buffer  # how many pixels to ignore along the x-axis
        self.y_buffer = y_buffer  # how many pixels to ignore along the y-axis
        self.patch_buffer = patch_buffer  # how many pixels to include before the tissue is detected
        self.patch_size = patch_size  # how much of a coordinate set to return
        self.patch_interval = 1 - patch_interval  # Fraction for how much the patches should overlap
        self.level = level  # smallest level of WSI from which to base your search on
        if self.wsi_path is None:
            raise RuntimeError('You must specify a WSI file path!')

    def get_targets(self):
        try:
            wsi = openslide.OpenSlide(self.wsi_path)
        except RuntimeError:
            raise RuntimeError('I was unable to read your WSI file: {}'.format(self.wsi_path))

        if self.level is not None:
            self.level = wsi.level_count - 1

        x, y = wsi.level_dimensions[self.level]
        img = wsi.read_region((0, 0), self.level, (x, y))
        # Change to color scale
        grey_img = img.convert('L')
        # Convert the image into numpy array
        np_grey = np.array(grey_img)
        # Identify where there is tissue
        # tuple where first element is rows, second element is columns
        idx = np.where(np_grey < self.threshold)
        patch_spacing = int(self.patch_size * self.patch_interval)
        previous_x_idx = 0

        for x_idx in range(idx[0].shape[0]):

            # Set the range so I stay in bounds of my slide
            min_x = max(idx[0][x_idx] - self.patch_buffer, 0)
            max_x = max(idx[0][x_idx] + self.patch_buffer, x)
            x_range = range(min_x, max_x)

            min_y = max(idx[1][x_idx] - self.patch_buffer, 0)
            max_y = max(idx[1][x_idx] + self.patch_buffer, 0)
            y_range = range(min_y, max_y)

            # If I am in between my target x and y's and I'm past my last point
            if idx[0][x_idx] in x_range and idx[1][x_idx] in y_range and idx[0][x_idx] > previous_x_idx:
                # Extract level 0 coordinates
                x_start = int(idx[0][x_idx] * wsi.level_downsamples[self.level])
                y_start = int(idx[1][x_idx] * wsi.level_downsamples[self.level])

                # Since the low magnification images are so much smaller, we have to loop through \
                # to get all of the high magnification image.
                x_end = x_start + wsi.level_downsamples[self.level] * self.patch_size
                y_end = y_start + wsi.level_downsamples[self.level] * self.patch_size

                # Make a copy of original y-start so I can loop over it for all the x
                y_start_org = y_start

                while x_start < x_end:
                    while y_start < y_end:
                        yield (x_start, y_start)
                        y_start = y_start + int(self.patch_size * self.patch_interval)
                    x_start = x_start + int(self.patch_size * self.patch_interval)
                    y_start = y_start_org
                previous_x_idx = idx[0][x_idx] + patch_spacing


if __name__ == '__main__':
    obj = WsiSample(wsi_path='163388.svs')
    counter = 0
    for i, j in obj.get_targets():
        counter += 1
    print('Found {} image patches'.format(counter))

