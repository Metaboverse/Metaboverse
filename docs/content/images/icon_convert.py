# Source: https://github.com/retifrav/python-scripts/blob/master/generate-iconset/generate-iconset.py
# Author: retifrav
# Available under MIT License
###
#MIT License
#
#Copyright (c) 2018 retif
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
###

# Instructions:
# Run as: `python icon_convert.py metaboverse_logo.png`

import os
import sys
import pathlib
import subprocess

if len(sys.argv) < 2:
    print("No path to original / hi-res icon provided")
    raise SystemExit

if len(sys.argv) > 2:
    print("Too many arguments")
    raise SystemExit

originalPicture = sys.argv[1]
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

# https://developer.apple.com/design/human-interface-guidelines/macos/icons-and-images/app-icon#app-icon-sizes
ListOfIconParameters = [
    IconParameters(32, 1),
    IconParameters(32, 2),
    IconParameters(64, 1),
    IconParameters(64, 2),
    IconParameters(256, 1),
    IconParameters(256, 2),
    IconParameters(512, 1),
    IconParameters(512, 2),
    IconParameters(1024, 1),
    IconParameters(1024, 2)
]

# generate iconset
for ip in ListOfIconParameters:
    subprocess.call(["sips", "-z", str(ip.width), str(ip.width), originalPicture, "--out", os.path.join(iconsetDir, ip.getIconName())])
    print(f"Generated {ip.getIconName()}")

# convert iconset to icns file
print('')
print('Converting to icns file...')
subprocess.call(["iconutil", "-c", "icns", iconsetDir, "-o", os.path.join(destDir, f"{fname}.icns")])
