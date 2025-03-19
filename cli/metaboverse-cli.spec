# -*- mode: python3.8 ; coding: utf-8 -*-
import os
import sys
import platform
sys.setrecursionlimit(5000)
block_cipher = None

output_name = 'metaboverse-cli-' + str(platform.system().lower())

# Add all required hidden imports
hidden_imports = [
    'scipy.special.cython_special',
    'scipy.special._cdflib',
    'scipy.spatial.transform._rotation_groups',
    'sklearn',
    'cmath',
    'chardet',
    'chardet.charset_normalizer',
    'numpy',
    'pandas',
    'networkx',
    'requests',
    'certifi',
    'pickle',
    'json',
    'zipfile',
    'pathlib',
    'importlib.util',
    'importlib.machinery',
    'importlib.abc'
]

# Add all required data files
datas = [
    (
        os.path.join(
            'metaboverse_cli',
            'source_url.txt'
        ),
        'metaboverse_cli'
    ),
    (
        os.path.join(
            'metaboverse_cli',
            'analyze',
            'data',
            'metabolite_mapping.pickle.zip'),
        os.path.join(
            'analyze',
            'data')
    ),
]

# Add all required binaries
binaries = []

# Add all required scripts
scripts = [
    os.path.join('metaboverse_cli', '__main__.py'),
    os.path.join('metaboverse_cli', '__init__.py'),
    os.path.join('metaboverse_cli', 'arguments.py'),
    os.path.join('metaboverse_cli', 'utils.py'),
]

a = Analysis(
    scripts,
    pathex=[
        'metaboverse_cli'
    ],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    runtime_hooks=[],
    excludes=[
        os.path.join('metaboverse_cli', 'test'),
        os.path.join('metaboverse_cli', '__test__.py'),
        os.path.join('metaboverse_cli', 'curate', 'test'),
        os.path.join('metaboverse_cli', 'curate', '__test__.py'),
        os.path.join('metaboverse_cli', 'mapper', 'test'),
        os.path.join('metaboverse_cli', 'mapper', '__test__.py'),
        os.path.join('metaboverse_cli', 'analyze', 'test'),
        os.path.join('metaboverse_cli', 'analyze', '__test__.py'),
        os.path.join('metaboverse_cli', 'analyze', 'test', 'HSA.zip'),
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False
)

pyz = PYZ(
    a.pure,
    a.zipped_data,
    cipher=block_cipher
)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name=output_name,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
