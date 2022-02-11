BoundingEllipsoidComputation 
============================================
This tool is for computing a bounding Ellipsoid around a set of points in n-dimensional space (for some value of n).

It bases heavily on GPL-v2ed code by Bojan Nikolic available from http://www.bnikolic.co.uk/blog/cpp-khachiyan-min-cov-ellipsoid.html, so this tool is available under GPLv2 as well.


Installation
============

Requirements
------------
- A moderately modern C++ and C compiler installed in a Unix-like environment, including the C++ library boost. Linux and MacOS should be fine.

To build the tool, run:

> cd src; qmake Tool.pro; make

The executable will be put into the src directory. Note that you may have to install some packages for the compilation process to succeed.

Input and output format
=======================

The input format is pretty self-explaining. Look at the example in the "examples" directory. Every point is on a separate line. 

The output of the tool is currently subject to change.

