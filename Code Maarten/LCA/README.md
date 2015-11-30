
## Welcome to the restyled LCA code of Maarten Vanhalst. ##


1. BUILDING
=============

This project now uses CMake to build and install. Perform the following steps:

	- $> cd build;
	- $> cmake ..
	- $> make (-j[numOfProcessors])
	- $> make install

[numOfProcessors] is the number of processors you want to use for building, this flag is optional.


If you have doxygen installed you can generate html documentation by doing

	- $> cd build;
	- $> make doc;

The documentation will be generated in the doc/ directory.

2. CODE CONVENTIONS
=====================
For the indentation we adopt the Kernighan & Ritchie style (K&R style) with an indentation with of 4 spaces. With the program "astyle" you can do this with
    - $> astyle --style=kr --indent=spaces=4 [yourfile.ext]

