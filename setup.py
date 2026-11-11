##########
# setup.py
##########

# Builds the Python module QuartetPatternCounterModule:  pip install .
# pybind11 is needed at build time; pip installs it from pyproject.toml.

import re
import sys
from pathlib import Path

from setuptools import setup, Extension
import pybind11

HERE = Path(__file__).resolve().parent

# The version is set only in QuartetPattern-counter-core.hpp.
version = re.search(r'#define\s+QPC_VERSION\s+"([^"]+)"',
                    (HERE / 'QuartetPattern-counter-core.hpp').read_text()).group(1)

if sys.platform == 'win32':
    compile_args = ['/std:c++17', '/O2']
    link_args = []
else:
    # No -march=native: the module must also run on CPUs other than the one that built it
    # (e.g. conda packages and wheels).
    compile_args = ['-std=c++17', '-O3', '-pthread']
    link_args = ['-pthread']

ext_modules = [
    Extension(
        'QuartetPatternCounterModule',
        ['QP_python_bindings.cpp', 'QuartetPattern-counter-core.cpp'],  # Both files!
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=compile_args,
        extra_link_args=link_args,
    ),
]

setup(
    name='QuartetPatternCounterModule',
    version=version,
    description='Site pattern frequencies of all quartets of sequences in a multiple sequence alignment',
    long_description=(HERE / 'README.md').read_text(encoding='utf-8'),
    long_description_content_type='text/markdown',
    author='Christoph Mayer',
    url='https://github.com/cmayer/QuartetPatternCounter',
    license='BSD-3-Clause',
    python_requires='>=3.8',
    install_requires=['numpy'],
    ext_modules=ext_modules,
    # Only the compiled module is installed. The .py files in this folder are tests and examples.
    py_modules=[],
    packages=[],
    zip_safe=False,
)
