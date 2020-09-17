# from https://github.com/retifrav/python-scripts/blob/master/generate-iconset/generate-iconset.py

import os
import sys
import pathlib
import subprocess

if len(sys.argv) < 3:
    print("No path to original / hi-res icon provided")
    raise SystemExit

if len(sys.argv) > 3:
    print("Too many arguments")
    raise SystemExit

originalPicture = sys.argv[1]
ouputDir = sys.argv[2]
if not (os.path.isfile(originalPicture)):
    print(f"There is no such file: {sys.argv[1]}")
    raise SystemExit

fname = pathlib.Path(originalPicture).stem
ext = pathlib.Path(originalPicture).suffix
destDir = pathlib.Path(originalPicture).parent

iconsetDir = os.path.join(destDir, f"{fname}.iconset")
if not (os.path.exists(iconsetDir)):
    pathlib.Path(iconsetDir).mkdir(parents=False, exist_ok=True)

class IconParameters():
    width = 0
    scale = 1
    def __init__(self,width,scale):
        self.width = width
        self.scale = scale
    def getIconName(self):
        if self.scale != 1:
            return f"icon_{self.width}x{self.width}{ext}"
        else:
            return f"icon_{self.width//2}x{self.width//2}@2x{ext}"

from PIL import Image
filename = r'logo.png'
img = Image.open(filename)
icon_sizes = [(16,16), (32, 32), (48, 48), (64,64), (128,128), (255,255)]
img.save('logo.ico', sizes=icon_sizes)
