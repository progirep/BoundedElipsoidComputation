BoundingEllipsoidComputation 
============================================
This tool is for computing a bounding Ellipsoid around a set of points in n-dimensional space (for some value of n).

It bases heavily on GPL-v2ed code by Bojan Nikolic available from http://www.bnikolic.co.uk/blog/cpp-khachiyan-min-cov-ellipsoid.html, so this tool is available under GPLv2 as well.

The tool executes Khachiyan's algorithm to compute the ellipsoid, followed by a straight-forward ellipsoid extension to ensure that all points are really contained in the ellipsoid.

The tool is an academic prototype by Ruediger Ehlers.

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

The output of the tool is currently unstable but self-explanatory.


Invoking the tool
=================

The tool takes a file name as parameter and some options. If the file name is not given, the tool will read from stdin instead.

The options are as followed:

* `--maxPrecision` followed by a double - Specifies the maximum precision (epsilon value) after which when the algorithm reaches it, the iterations in Khachiyan's algorithm are stopped. Set this to 0.0 to abort when the algorithm's iterations do not change the elipsoid any more. The default value is 0.00001
* `--maxIterations` followed by an integer - Specifies how many iterations of Khachiyan's algorithm are executed at most. This option and the preciding one both specify termination criteria - the execution of Khachiyan's algorithm continues while none of the criteria indicate termination.
* `--printVolumeEvery` followed by an integer. By default, the tool is quite during the main computation. If this option is given with an integer denoting a number of iterations of the algorithm's main loop, every so many iterations, the volume of the current Ellipsoid is printed. Note that this is actually expensive.
* `--skipfirstNLines` followed by an integer - Ignores the first specified number of lines from the input. This may be useful for feeding the tool with the output of the tool `qhull`, which will put a header to the points file. This is currently untested.
* `--doNotInstallSIGINTHandler` - By default, the tool will finish the current iteration of Khachiyan's algorithm and compute and display the ellipsoid computed so far then the tool receives a SIGTERM signal. If this is not desired and the tool should rather terminate immediately without providing a result when receiving this signal, this option can be used.







