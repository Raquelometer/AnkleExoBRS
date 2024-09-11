function [g, data, tau2] = computeAnkleExoBRSwithExo()
%% Function to compute baseline backward reachable sets (without exoskeleton)

%% Load files
model = 'AnkleExoBRS/Models/OM_model.mat'
problemParams = 'AnkleExoBRS/Params_OM_gcExo_100_pctMT_100_pctRTD.mat'

load(model)
load(problemParams)

%% Construct dynSys object
x0 = [pi/2, 0, 0];
dPend = LiftedInvPendExo(x0, model, problemParams)

%% grid
N = [225; 225; 225];         % Number of grid points per dimension

% Partial grid if forward velocity
% grid_min = [0.65, -0.02*dPend.vscale, dPend.torqueMin*dPend.tauscale]; % Lower corner of computation domain
% grid_max = [1.8; 2.85*dPend.vscale; dPend.torqueMax*dPend.tauscale];    % Upper corner of computation domain

% Partial grid for negative velocity
grid_min = [1.51, -2.85*dPend.vscale, dPend.torqueMin*dPend.tauscale]; % Lower corner of computation domain
grid_max = [2.56; 0.02*dPend.vscale; dPend.torqueMax*dPend.tauscale];    % Upper corner of computation domain

% TODO add exo grid

% Uncomment for full grid
% grid_min = [0.65, -2.85*dPend.vscale, dPend.torqueMin*dPend.tauscale]; % Lower corner of computation domain
% grid_max = [2.56; 2.85*dPend.vscale; dPend.torqueMax*dPend.tauscale];  

g = createGrid(grid_min, grid_max, N);

% dimensionless time vector
t0 = 0;
tMax = 1.4/dPend.omega;
dt = 0.7/dPend.omega;
tau = t0:dt:tMax;

%% Target sets

TS = problemParams.targetSets;

% TODO can't assume exo is included!!!, add scaling condition
data0 = shapeRectangleByCorners(g,[TS.ankleEq(1)- TS.thetaAnkRange; -0.01*dPend.vscale; -1*dPend.tauscale],[TS.ankleEq(1) + TS.thetaAnkRange; 0.01*dPend.vscale; 1*dPend.tauscale]);
data1 = shapeRectangleByCorners(g,[TS.toeEq(1) - TS.thetaToeRange; -0.01*dPend.vscale; (TS.toeEq(3) - TS.torqueToeRange)*dPend.tauscale],[TS.toeEq(1) + TS.thetaToeRange; 0.01*dPend.vscale; (TS.toeEq(3) + TS.torqueToeRange)*dPend.tauscale]);

data0 = shapeUnion(data0,data1);

% Constraints
% constraintMat = generateLiftedInvPendExoCoPConstraint(g.vs, g.shape, dPend);
CoPConstraintMat = generateLiftedInvPendExoCoPConstraint(g.vs, g.shape, dPend);
GRFConstraintMat = generateLiftedInvPendExoGRFConstraint(g.vs, g.shape, dPend);
constraintMat = min(CoPConstraintMat,GRFConstraintMat);

% Friction
mu = 1; % coefficient of static friction
FrictionConstraintMat = generateLiftedInvPendFrictionConstraint(g.vs, g.shape, dPend, mu);
constraintMat = min(constraintMat,FrictionConstraintMat);
%% HJB Solver Params

% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = dPend;
schemeData.accuracy = 'high'; %set accuracy
schemeData.uMode = 'min';

% Compute value function
HJIextraArgs.obstacleFunction = constraintMat;
%HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 1;
HJIextraArgs.visualize.plotColorVS0 = 'blue';
HJIextraArgs.visualize.plotAlphaOS = 0.3;
HJIextraArgs.visualize.plotColorOS = 'green';
HJIextraArgs.visualize.plotColorOF = 'green';
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = false; %delete previous plot as you update
% HJIextraArgs.compMethod = 'minVOverTime';
HJIextraArgs.compMethod = 'minVWithV0';

HJIextraArgs.quiet = 1;

[data, tau2, ~] = ...
  HJIPDE_solve(data0, tau, schemeData, 'none', HJIextraArgs);
