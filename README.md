# Ising-2D
Fortran code to simulate the Ising model in 2 dimensions. 
This code was developed for the [Cooperative and critical phenomena](https://ifisc.uib-csic.es/master/programme-syllabus/cooperative-and-critical-phenomena/) subject imparted in the [Complex Systems master course](https://ifisc.uib-csic.es/master/) from the [IFISC](https://ifisc.uib-csic.es/es/).
The simulation includes the computations for the magnetization, energy and correlation of the system. A final report summarizing the dynamics of the system and the main outputs is provided in the **simulation-ising-model.pdf** file.

## Ising model
The Ising model was proposed by Wilhelm Lenz and Ernst Ising in 1920 to study the behavior of ferromagnetic materials, and today it has become one of the major paradigms of statistical mechanics. 
In 1944 Lars Onsager managed to find the analytical solution to the problem, but the development is arduous and cumbersome. 
Therefore, in this work we will analyze the behavior of this system in its 2-dimensional version through simulations. 
The various magnitudes obtained will have associated errors of different nature, for example due to the finitude of the systems studied.


## Code files
- dranxor.f90. Random number generator in Fortran90 from "Generation of Gaussian distributed random numbers by using a numerical inversion method", R. Toral, A. Chakrabarti. Computer Physics Communications, 74 (1993) 327-334. Please give credit to the authors of the code R. Toral, A. Chakrabarti, and properly cite them when using or disseminating it. This repository includes it under the free use/dissemination statement of the authors. 
- isingtest5.f. MC algorithm using the Metropolis and Wolff methods for the 2D Ising model, to compute both magnetization and energy
- isingtest8.f. MC algorithm using the Metropolis and Wolff methods for the 2D Ising model, to compute correlation
- parameters. Parameters used for the simulations
- plots.ipynb. Very brief Ipython code to generate plots from the output files


## Contributors
Patrick Sánchez Galea ([Kaggle profile](https://www.kaggle.com/saga21))
