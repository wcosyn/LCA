
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


2. FOLDER INFO
================

    - ./src/ contains the source files from and the main (mainop) left behind by Maarten. This is a large main file for all different calculations.
      It could be more convenient to split this up into executables calculating just one thing, with a suiting filename.

    - ./main/ contains smaller source files containing a main. The philosophy is to make small main files dedicated to calculating one thing.
      (ob momentum distr., tb rel. momentum distr, tb c.m. momentum distr. , ...)

    - ./include/ contains the header files

    - ./data/ contains the data that is used by the executables, such as saved moshinsky brakets in ./data/mosh

    - ./manual/ contains the manual pdf

    - ./results/ directory for the calculated results. The new results are in ./results/new_results2017/

    - ./scripts/ contains scripts Maarten left behind. I doubt that there will be scripts that will be of much use/compatible with latest version of the code.
   
    - ./test/ directory with source file containing a main() written for testing parts of the code

    


3. CODE CONVENTIONS
=====================
For the indentation we adopt the Kernighan & Ritchie style (K&R style) with an indentation with of 4 spaces. With the program "astyle" you can do this with
    - $> astyle --style=kr --indent=spaces=4 [yourfile.ext]


