# Path-Integral-Monte-Carlo
This repository contains codes for solving Path-Integral by Monte-Carlo methods. The codes are written in both Python and Julia. Due to faster computation times and ability to do multithreading Julia is recommended. 

The following are the required packages for Python:
* Numpy
* Matplotlib
* random
* time

Also, these are written with Python 3.8.5. Thus it is recommended to use this or higher version to run.

The Julia codes are written with Julia version 1.5.0. The following packages in Julia are used:
* Plots
* PyPlot
* Random
* LaTeXStrings

These codes use multithreading to do computations quickly, thus use the following command to increase the number of threads used:  
`export JULIA_NUM_THREADS=4`  
for 4 threads. This command is for bash shell, and remains active till the terminal window is closed.  
To check the number of threads in Julia, one can print:  
`Threads.nthreads()`  

There are two folders one for Harmonic oscillators and another for Anharmonic oscillators. The folders contains codes for each of the various parameters to compute.
They can be used to make the plots shown in the Computational_physics_report.pdf file.


