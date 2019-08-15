City-Scale CFD
===
This repository contains several packages of our CFD model implementations used in a few recent papers.

These models can be used to simulate the city-scale thermal and wind environment in very fine resolution, such as urban heat island circulation (UHIC) simulation, including both the quasi-steady UHIC and UHIC evaluation, which is an important problem in urban envirionment study. Some key features in these models are:
1)	Governing equations based on new coordinate.
The governing equations (including the momentum equation, energy equation and turbulence equations) are obtained with KRB coordinate transformation. This can be implemented by DEFINE_SOURCE macro.
2)	Absorbing layer.
An absorbing layer is added below the top boundary to avoid the effect of spurious waves.
3)	Mixed porous approach.
The city area is treated as porous medium.
4)	Daily cycle temperature boundary condition.
An analytical surface temperature formula was derived and implemented in these models after simplifying longwave radiation and sensible heat flux terms.

Usage
===
ANSYS FLUENT 14.0 is used for running these models.

Five packages are included. Each package includes a source code file (C function that can be either interpreted or compiled in ANSYS FLUENT) "\*.c", a journal file (a simple text file containing Fluent Text User Interface (TUI) commands) "\*.jou" and a mesh file "\*.msh". Details of these packages are:<br>
  2D: 2D quasi-steady urban heat island circulation (UHIC) simulation<br>
  2DPor: 2D quasi-steady UHIC in which the city is treated as porous medium<br>
  2D12hr: 2D UHIC evolution including the city effect<br>
  3D: 3D quasi-steady UHIC simulation<br>
  3DPor: 3D quasi-steady UHIC in which the city is treated as porous medium

To use one package, one may download it and load the journal file. One may want to change the mesh file if needed. The C functions are called from journal files.  

Citation
===
If you make use of these codes, please cite relevant publications from the following:

  [1] Wang, X., & Li, Y. (2016). Predicting urban heat island circulation using CFD. Building and Environment, 99, 82-97<br>
  [2] Wang, X., Li, Y., & Hang, J. (2017). A combined fully-resolved and porous approach for building cluster wind flows. Building Simulation, 10(1), 97-109<br>
  [3] Wang X., Li Y., Wang K., Chan PW. and Yang X. (2017), A simple daily cycle temperature boundary condition for ground surfaces in CFD predictions of urban wind flows, Journal of Applied Meteorology and Climatology, 56(11), 2963-2980.

Copyright
===
These packages are provided for use for academic purposes only. Any commercial use without authorization is not permitted and for commercial uses, please contact Yuguo Li at liyg@hku.hk.
