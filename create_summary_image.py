import os
import argparse
import sys
import pwd
import time
import subprocess
import re
import shutil
import glob	
import openslide
import numpy as np
from PIL import Image, ImageDraw
import math
import tensorflow as tf
import io
from shapely.geometry import Polygon, Point, MultiPoint
from shapely.geometry import geo
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
import tensorflow as tf
import io
#from dataset_utils import * 

svs='/projects/shart/digital_pathology/data/TCGA/test_svs/test/TCGA-A1-A0SM-01Z-00-DX1.AD503DBD-4D93-4476-B467-F091254FDF78.svs'
OSobj = openslide.OpenSlide(svs)
inlevel=3
level=inlevel
input_level=0
patch_sub_size_x=OSobj.level_dimensions[inlevel][0]
patch_sub_size_y=OSobj.level_dimensions[inlevel][1]

img_patch = OSobj.read_region((0,0), level, (patch_sub_size_x, patch_sub_size_y))
multi_factor=OSobj.level_downsamples

poly_included=[]
poly_excluded=[]
fobj = open("k.txt")
for i in fobj:
	i = i.strip()
	p = i.split("_")
	name=p[0]
	#name = name[:-1]
	x1=int(p[1])
	x2=x1+512
	y1=int(p[2])
	y2=y1+512

	if x1 == "0":
		x1 = int(x1)
	else:        
		x1=(int(x1)*multi_factor[input_level])/multi_factor[level]
		x1=int(x1)
	if y1 == "0":
		y1 = int(y1)
	else:    
		y1=(int(y1)*multi_factor[input_level])/multi_factor[level]
		y1=int(y1)

	x2=(int(x2)*multi_factor[input_level])/multi_factor[level]
	x2=int(x2)
	y2=(int(y2)*multi_factor[input_level])/multi_factor[level]
	y2=int(y2)
	#if 'included' in i:
	poly_included.append(Polygon([(x1,y1),(x2,y1),(x2,y2),(x1,y2),(x1,y1)]))
	#else:
	#poly_excluded.append(Polygon([(x1,y1),(x2,y1),(x2,y2),(x1,y2),(x1,y1)]))

f, ax = plt.subplots(frameon=False)
f.tight_layout(pad=0, h_pad=0, w_pad=0)
ax.set_xlim(0, patch_sub_size_x)
ax.set_ylim(patch_sub_size_y, 0)
ax.imshow(img_patch)
#for j in range(0,len(poly_excluded)):
	#patch1 = PolygonPatch(poly_excluded[j], facecolor=[0,0,0], edgecolor="red", alpha=0.4, zorder=2)
	#ax.add_patch(patch1)
for j in range(0,len(poly_included)):
	patch1 = PolygonPatch(poly_included[j], facecolor=[0,0,0], edgecolor="red", alpha=0.3, zorder=2)
	ax.add_patch(patch1)        
ax.set_axis_off()
DPI = f.get_dpi()
plt.subplots_adjust(left=0, bottom=0, right=1, top=1,wspace=0, hspace=0)
f.set_size_inches(patch_sub_size_x / DPI, patch_sub_size_y / DPI)
f.savefig("sample.png", pad_inches='tight')
fobj.close()
