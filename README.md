# AnkleExoBRS
All dependencies are in MATLAB
## Target sets code requires:
  1. [SeDuMi](https://sedumi.ie.lehigh.edu/?page_id=58)
 
  2. [N-dimensional Convex Polyhedra 
toolbox](https://www.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra)

## BRS computation requires
   1. [Hamilton-Jacobi Reachability Toolbox](https://github.com/HJReachability/helperOC)

## Computing stabilizable regions with precomputed target sets
Precomputed target sets for all of the scenarios in the paper are saved in the TargetSets directory. Models are in the Models directory.
   1. Download and install the HJB toolbox
   2. Move the Dynamics/LiftedInvPendExo directory into the dynsys directory that comes with the toolbox
   3. Make sure the functions in Dynamics/Constraints and in Models/ModelUtils are added to your path.
   4. The computeNominalBRS.mat and computeAnkleExoBRS.mat scripts will compute the BRS. Select the model you want (i.e. young male = YM_model.mat) and the desired problem parameters. 
