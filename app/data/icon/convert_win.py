# from https://github.com/retifrav/python-scripts/blob/master/generate-iconset/generate-iconset.py

import os
import sys
import pathlib
import subprocess
import shutil
from PIL import Image

if len(sys.argv) < 3:
    print("No path to original / hi-res icon provided")
    raise SystemExit

if len(sys.argv) > 3:
    print("Too many arguments")
    raise SystemExit

originalPicture = sys.argv[1]
ouputDir = sys.argv[2]

originalPicture = "/Users/jordan/Desktop/projects/Metaboverse/app/data/icon/png/icon_512x512.png"
ouputDir = "/Users/jordan/Desktop/projects/Metaboverse/app/data/icon/win"

if not (os.path.isfile(originalPicture)):
    print(f"There is no such file: {sys.argv[1]}")
    raise SystemExit

fname = pathlib.Path(originalPicture).stem
ext = pathlib.Path(originalPicture).suffix
destDir = pathlib.Path(ouputDir)

iconsetDir = os.path.join(destDir, "metaboverse_logo.ico")

if not (os.path.exists(iconsetDir)):
    pathlib.Path(iconsetDir).mkdir(parents=False, exist_ok=True)

img = Image.open(originalPicture)
icon_sizes = [(16,16), (32, 32), (48, 48), (64,64), (128,128), (255,255)]
img.save("metaboverse_logo.ico", sizes=icon_sizes)
shutil.move('metaboverse_logo.ico', 'win/metaboverse_logo.ico')
