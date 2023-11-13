# KEMPO2: Space Plasma Simulator

*KEMPO2 (Kyoto university ElectroMagnetic Particle cOde in 2D)* is a two-dimensional Particle-in-Cell code with Graphical User Interface for the analysis of whistler-mode instability. The simulation interface is designed with App Designer in MATLAB environment. This simulator was developed by Kenta Hirahra, as the final project of the Bachelor's studie.

# Requirements

MATLAB 2019a or higher.

# Installation

```
$ git clone https://github.com/kenta-hirahara/kempo2_sotsuron.git

```

# Code execution

Go to the folder in MATLAB and run 
```
$ KEMPO2

```
in the command window. Then the UI shows up and you can set parameters as you wish.

# Note

 - For meaningful results, it is imperative that the grid size be a minimum of 64 squared units in the xy plane, accompanied by a particle count of 128. This stipulation arises from the application of a Gaussian distribution to the particles, requiring a sufficient grid resolution for accurate simulation outcomes.
 - Calculation resouces for velocity distribution, E\*J and each dispersion relation are quite large. So please remove the ticks from the check boxes in the Panel Options if they're unnecessary.
 - We recommend ntime to be a power of two for Fourier Transformations.
