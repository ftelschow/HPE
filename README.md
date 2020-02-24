# HPE
This repository implements the bootstrapped Hermite Projection estimator
and other estimators for Lipschitz Killing curvatures of random fields.
Moreover, it reproduces the simulation results and data analysis from
https://arxiv.org/abs/1908.02493.

## Table of contents
* [Introduction](#introduction)
* [Set Up](#setup)
* [Folder Structure](#folderstruct)
    * [3rdparty](#3rdparty)
    * [code](#MatlabFunctions)
      * [EEC](#EEC)
      * [LKC](#LKC)
      * [RandomFields](#RandomFieldGeneration)
      * [auxilaryFcn](#Auxfunctions)
    * [scripts](#scripts)
 

## Introduction <a name="introduction"></a>
This git-repository implements and aims to reproduce the
simulation results and data analysis from "Estimation of Expected Euler
Characteristic Curves of Nonstationary Smooth Gaussian Random Fields" (Telschow et al)
using Matlab. Moreover, it will also host simulations for an upcoming NeuoImage article
analyzing the LKC estimators in 3D.mex EulerCharCrit_c.cpp

These estimators will be soon included into the RFTtoolbox https://github.com/sjdavenport/RFTtoolbox,
such that the latter Toolbox should be used for future research, since we will keep only the RFTtoolbox
fully up-to-date in the future and new developements will appear there. Please cite that Toolbox appropriately,
in case you are using our work.

If you find bugs or you have problems setting up the repository, feel free to get in touch with me ftelschow(AT)ucsd.edu.

## Set Up <a name="setup"></a>
This repository requires spm12, which is an open source
neuroimaging matlab toolbox available under https://www.fil.ion.ucl.ac.uk/spm/software/spm12/.
Please download it and follow the corresponding install instructions before proceeding.

If you not want to run any of the generate_profile.m files from the scripts folder, which we
recommend and will set up after changing the appropritate paths in the script, go to the folder
**/HPE/code/EEC/csource** and run the command
```
mex EulerCharCrit_c.cpp
```
in order to compile the C++ code.


## Folder Structure <a name="folderstruct"></a>
### 3rdparty <a name="3rdparty"></a>
This folder collects functions written by 3rd parties. Mostly they are used to
transform matlab arrays into .nii files to use spm12.
### code <a name="code"></a>
This folder contains our implemented functions and has a substructure according
to the purpose of the function, which is detailed below.
#### EEC <a name="EEC"></a>
This folder contains functions to estimate the expected euler characteristic curves
of random fields. It uses a critical value approach and partially is coded in C++ in
order to speed up the code.
#### LKC <a name="LKC"></a>
This folder implements different estimatorsmex EulerCharCrit_c.cpp for the LKCs from the literature and the
new Hermite Projection Estimators (HPE).
#### RandomFields <a name="RandomFieldGeneration"></a>
This folder provides different functions to generate random fields over 1-,2- and 3-D domains.
#### auxilaryFcn <a name="Auxfunctions"></a>
This folder contains small functions helping with the data analysis.

### scripts <a name="scripts"></a>
This folder contains subfolders for the two different articles dealing with LKC estimation.
Each folder contains a setup script called generate_profile.m, which needs to be run first before
any of the other scripts in the folder. It will automatically produce the neccessary subfolders and
compile the C++ code. **Don't forget to change the paths in this script according to your folder structure!**
Afterwards all other scripts, which reproduce the simulation results, data analysis and
figures from the paper can be run.

README.md
