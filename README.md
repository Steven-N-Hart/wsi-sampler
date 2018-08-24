This class will scan the highest level of a whole slide image to identify regions with sufficient tissue


Args:
```
        wsi_path: '/path/to/wsi.svs' [Required]
        threshold = Value at which one considers empty space [200]
        x_buffer = How many pixels to ignore along the x-axis [50]
        y_buffer = How many pixels to ignore along the y-axis [50]
        patch_buffer = How many pixels to include before the tissue is detected [5]
        patch_size = How large of a coordinate set to return [300]
        patch_interval = Fraction for how much the patches should overlap [0.1]
```

Usage:
```
        import WsiSample
        import openslide

        wsi_path='/path/2/wsi.svs'
        patch_size = 300
        level = 0

        obj = WsiSample(wsi_path='/path/2/wsi.svs', patch_size=patch_size)
        OSobj = openslide.OpenSlide(wsi_path)

        for x,y in obj.get_targets(): # x and y are coordinates to extract from your OpenSlide Object
            img_patch = OSobj.read_region((x,y), level, (patch_size, patch_size))
 ```