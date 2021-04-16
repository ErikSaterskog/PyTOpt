# PyTOpt
This repository contains all of the files necessary for execution of the topology optimisation program written by Daniel Pettersson and Erik Säterskog.

## Installation
  The program is run by executing one of the files in the "Test_Examples" folder or creating your own structure in the same manner.

  Packages required:

  scipy\
  numpy\
  visvis\
  pyvtk\
  calfem-python\
  sphinx==1.6.6

  For the Meshgenerator to work must GMSH.exe be in the same folder as the test file. GMSH.exe can be downloaded from the Repository or at https://gmsh.info/.

  To install the package write: 'pip install git+https://github.com/ErikSaterskog/Thesis_Repo.git --user' in your terminal.\ 
  Be sure that setuptools is installed before installing this package.\
  DISCLAIMER: Git may must be installed and added to PATH. In the future will the package hopefully be uploaded to pypi.org and thereafter will this not be an issue.

## Lincense
  Copyright (c) 2021 Daniel Pettersson\
  Copyright (c) 2021 Erik Säterskog
  
  PyTOpt is under the terms of the GNU General Public License as published by the Free Software Foundation. This means that it is free to distribute and modify. However, this    program is without warranty. For more information about the terms and conditions read the LICENSE file or go to <http://www.gnu.org/licenses/>.

  
## References

  [Svanberg, K. (1987). The Method of Moving Asymptotes – A new method for structural optimization. International Journal 
  for Numerical Methods in Engineering 24, 359-373. doi:10.1002/nme.1620240207](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620240207)

  Svanberg, K. (n.d.). MMA and GCMMA – two methods for nonlinear optimization. https://people.kth.se/~krille/mmagcmma.pdf 

  See master's thesis [ref] for more references.
