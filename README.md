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
   2. Move the `Dynamics/LiftedInvPendExo` directory into the dynsys directory that comes with the toolbox
   3. Make sure the functions in `Dynamics/Constraints`, `Dynamics/DynUtils`, and `Models/ModelUtils` are added to your path.
   4. The `computeNominalBRS.mat` and `computeAnkleExoBRS.mat `scripts will compute the BRS. 
 Within the scripts, load the files corresponding to the model you want and the desired problem parameters. 
  
For example for a Young Male with full strength, load `/Models/YM_model.mat` and problem parameter file `TargetSets/computedSets/Params_YM_noExo_100_pctMT_100_pctRTD.mat`.

You should also make sure to uncomment the grid portion for the forward or backward velocity grid, depending on which portion of the BRS you want to compute.
