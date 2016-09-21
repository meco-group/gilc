# gILC
gILC was originally developed by Marnix Volckaert in 2012. It is currently maintained and adapted for new research topics by Armin Steinhauser within the [MECO research team](https://www.mech.kuleuven.be/en/pma/research/meco/).

Make sure to get [CasADi](http://casadi.org/) (at least v3.0) before using gILC.

# Notes
The originally published software (v1.3) can be found [here](https://set.kuleuven.be/optec/Software/gilc-generic-iterative-learning-control-for-nonlinear-systems). In case this website is not accessible, the content is packed in the *mirror-v1.3* folder.

# Version history
- v1.4, 09/2016
  - Changed to CasADi v3.0 syntax
  - Expecting separate CasADi installation and abandoning the install-paradigm
  - Move to Github repo


- v1.3, 05/2012
  - Updated manual document


- v1.2
  - Included option to provide the number of additional states (which won't be set equal to the initial state value)
  - Changed the default weight on tracking error to support multiple outputs
  - Fixed a bug that misinterpreted multiple references (for MIMO systems)
  - Removed additional constraints for ending the motion at rest using the relative degree
  - Changed the default weight on tracking error to 1e8
  - Updated tutorials


- v1.1
  - Updated tutorials 1,2 and 3


- v1.0
  - Initial version
