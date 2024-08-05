# Turbulence anisotropy in a streamline coordinate system
This directory contains codes for computing geometric properties of the anisotropic Reynolds stress tensor.  
Eigenvalues are computed and used to define the barycentric map (Banerjee et al., [2007](https://www.tandfonline.com/doi/full/10.1080/14685240701506896)), which is a variant of the Lumley triangle.
Kite-shaped regions at the edges of the barycentric map (Stiperski and Calaf, [2018](https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.3224)) are defined to select Reynolds stress tensors corresponding to the three limiting states of turbulence: 
one-component, two-component axisymmetric and isotropic turbulence.
Eigenvectors are computed and their directions analysed in relation to the streamline coordinates defined in the direction of the mean wind vector.

The corresponding publication is 'Interpreting turbulence anisotropy in a streamline coordinate system' (just submitted to JGR: Atmosphere).
In the directory 'Visualization', the reader finds the codes for reproducing the figures of the paper.
