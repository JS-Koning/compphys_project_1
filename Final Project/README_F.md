# Project 1: Molecular dynamics simulation of Argon atoms
last updated: 27/02/2021
<br>

This is the first project for Computational Physics (AP3082) year 2020/2021, T-course of the master Applied Physics at the Delft University of Technology.
The code is written for Python 3.8.6 using JupyterLab v2.2.9 web-based userinterface. 
<br>

This program plots the relative distance, force and energy of N argon particles with a maximum 3 dimensions.
The potential of the argon atoms is based on the Lennard-Jones potential and thus, so is the force. 
<hr>
To run the program: <br>
optional: edit parameters <br>
run skeleton.py (for supported versions of Python please see <strong>Installation</strong> 

<hr> 
<br>
<br>

## functionality


   - Dimensional, and non-dimensional SI units
       - Currently, all inputs and outputs are in either dimensional or non-dimensional form, no combination is supported
   - Periodic or non-periodic boundary conditions
   - Initial random locations
       - Only inside predefined box size
   - Initial 'random' velocities 
       - loosely based on 3D equipartiton theorem, changing dimensions in parameters does not change random velocity algorithm
   - Plots of total distance between atoms
   - Plots of total energy
   - Plots of kinetic energy
   - Plots of potential energy
   - Choice of Euler method numerical approximation or Verlet algorithm
   - Potential (and thus force) based on Lennard-Jones potential


## Installation
Usage of JupyterLab v2.2.9 is recommended, different versions might have reduced functionability/performance. For an installation manual please visit the [Jupyterlab documentation](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html).
Since this code is written using Python 3.8.6, it is likely that most features will work for various IDE's running Python 3.8.6.

## Acknowledgements
The following persons are the instructors of the course and supplied theoretical material:
<ul>
<li>Dr. Michael Wimmer</li>
<li>Dr. Anton Akhmerov</li>
</ul>
The following persons are the technical assistants of this course and are available for intermittent grading, feedback and questions: 
<ul>
<li>Helene Spring</li>
<li>Marta Pita Vidal</li>
<li>André Melo</li>
<li>Anastasiia Varentcoval</li>
<li>David van Driel</li>
</ul>
A public forum is available during this course where anyone related to this course can ask questions, and provide insights.</br>
JupyterHub has been supplied by the Delft University of Technology to allow cloud computing, without the installment of any software.

## Authors
- Abi Kanagaratnam
- Jim Koning

