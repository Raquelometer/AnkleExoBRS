function [TS_array, theta_lb, theta_ub] = generateExoTargetSets(model,modelLabel,alphaRTD, alphaMT, muGrav, motor_lim)
%"generateExoTargetSets" Function to compute box-shaped target sets
%
%   Function takes in a model that includes mass, height, and bounds on RTD
%   and torque. The user should supply the desired scaling factors alphaRTD and
%   alphaMT and the percentage of gravity compenstion muGrav, and the exoskeleton
%   motor limit.
%
%   Function output is an 8x2 array of box coordinates centered at the
%   ankle and toe equilibrium points.
plotting = true;
saveData = true;
load(model);
m = model.m;
h = model.h;
g = 9.8;
uMin = model.uMin*alphaRTD;
uMax = model.uMax*alphaRTD;

tau_min = model.tauMin*alphaMT;
tau_max = model.tauMax*alphaMT;


[l, lf, mf, ank, m, ~, ~] = proportionallyEstimatedParams(m, h);

%% Exo params
if muGrav ~= 0
    Kp = 0;
    Kd = 0;
    exoLabel = 'gcExo'
else
    Kp = m*g*l;
    Kd = 0.3*sqrt(m*l*l*Kp);
    exoLabel = 'PDExo';
end


%% Compute gravity compensation boundaries
theta_min = acos(ank/l);
theta_max = pi - acos((lf-ank)/l);


heel_torque = m*9.8*l*cos(theta_min)

toe_torque = m*9.8*l*cos(theta_max)

%% Plot relevant boundaries as sanity check

theta_vals = linspace(0, pi, 1000);

Kx = max(min((muGrav*m*g*l*cos(theta_vals) + Kp*((pi/2)-theta_vals)), motor_lim),-motor_lim);
Tau_vals = m*g*l*cos(theta_vals) - Kx;

theta_lb = theta_min - (pi/2);
theta_ub = theta_max - (pi/2);

theta_lb_exo_sat = acos(motor_lim/(m*g*l)) - (pi/2);
theta_ub_exo_sat = acos(-motor_lim/(m*g*l)) - (pi/2);

% theta_middle_exo_sat = theta_lb + (theta_ub_exo_sat - theta_lb)/2;

figure(1)
hold on
plot(theta_vals - (pi/2), Tau_vals);
yline(tau_min, 'Label', 'Torque Lower Bound')
yline(tau_max, 'Label', 'Torque Upper Bound')
xline(theta_min - (pi/2), 'Color','blue','LineWidth',1.5, 'Label', 'Heel','LabelVerticalAlignment','middle')
xline(theta_max - (pi/2), 'Color','blue','LineWidth',1.5, 'Label', 'Toe','LabelVerticalAlignment','middle')
% xline(theta_lb,'Color','green','LineWidth',1.5)
% xline(theta_ub, 'Color','green','LineWidth',1.5)
xline(theta_ub_exo_sat, 'Color','magenta','LineWidth',1.5,'Label','Exo saturation point', 'LabelVerticalAlignment','middle')
xline(theta_lb_exo_sat, 'Color','magenta','LineWidth',1.5,'Label','Exo saturation point', 'LabelVerticalAlignment','middle')
xlim([-0.2 .3]);
ylim([-150 50]);
hold off

% %% Adjust theta_lb and ub if there is enough torque
% theta_lb = max(theta_lb, theta_min - (pi/2));
% theta_ub = min(theta_ub, theta_max - (pi/2));
%% Set up 3d Plots
if plotting

    x1 = linspace(-1,1, 320);
    x2 = linspace(-1,1,320);
    x3 = linspace(-70,50,600);
    
    [X1, X2, X3] = meshgrid(x1, x2, x3);
end
%% Compute CLF at ankle 
eq_angle = 0;
current_eq = [eq_angle, 0, -m*g*l*sin(eq_angle)]';

% Linearized system without exo
n = 3;
A = [0, 1, 0; 0, -Kd/(m*l*l), 1/(m*l*l); 0, 0, 0];
B = [0;0;1];


% SDP
Q = sdpvar(n,n); 

Y = sdpvar(1,n);


eps = 1; 

Consts = [[Q>=eps*eye(n)]; [Y'*B'+Q*A' + A*Q + B*Y]<=0*eye(n)];

% Obj = ([trace(Q)]); % multiply Q with a scaling diagonal
Obj = ([]); % multiply Q with a scaling diagonal


ops = sdpsettings('solver', 'sedumi');

diagnostics = solvesdp(Consts, Obj, ops);

if diagnostics.problem == 0

    disp('Feasible: CLF found')
    Qval = double(Q); %use "double" to extract solution
    Yval = double(Y);
    P = inv(Qval);
    K = Yval*P

end
R = chol(P);

%% Find level
sway_velocity = 0.01;

% Construct constraint polytope 
A_con = [K; -K; 1 0 0; -1 0 0; 0 1 0;0 -1 0;0 0 1;0 0 -1];
b_con = [uMax; -uMin; theta_ub; -theta_lb; sway_velocity; sway_velocity; tau_max; -tau_min;]
verts = lcon2vert(A_con,b_con);

% Select level set candidate (max ellipsoid contained in the constraints)
level_candidate = computeLevelCandidate(A_con, b_con, Qval, R, current_eq)


% Compute the minimum circumscribed box 
[axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate)

% convert to vertices
ellipseBoxVerts = boxVerts3d(axes_max, current_eq);
% figure(3)
% hold on 
% scatter3(verts(:,1), verts(:,2),verts(:,3), 50, 'red', 'filled');
% 
% scatter3(ellipseBoxVerts(:,1), ellipseBoxVerts(:,2),ellipseBoxVerts(:,3), 50, 'green', 'filled');
% 
% hold off

% convert vertices to linear constraints using vert2lcon
[A_ts, b_ts] = vert2lcon(ellipseBoxVerts);


% Compute the vertices of the constraint polytope intersected with the
% bounding box of the max inscribed elliposid
verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts]);

target_tol = 0.0001;

while size(verticesIntersection,1) ~= 8  

    level_candidate = level_candidate*0.95;
    [axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate);
    ellipseBoxVerts = boxVerts3d(axes_max);
    [A_ts, b_ts] = vert2lcon(ellipseBoxVerts);
    verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts]);

end


while max(abs(verticesIntersection - ellipseBoxVerts),[], 'all') > target_tol 
    level_candidate = level_candidate*0.95;
    [axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate);
    ellipseBoxVerts = boxVerts3d(axes_max,current_eq);
    [A_ts, b_ts] = vert2lcon(ellipseBoxVerts);
    verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts]);
end


TS_array = ellipseBoxVerts;
if plotting
f = @(x1, x2, x3, S) S(1,1).*x1.*x1 + S(1,2).*x1.*x2 + S(1,3).*x1.*x3 + S(2,1).*x1.*x2 + S(2,2).*x2.*x2 + S(2,3).*(x3.*x2)...
    + S(3,3).*x3.*x3  + S(3,1).*x1.*x3 + S(3,2).*x2.*x3;
V = f(X1 - current_eq(1), X2 - current_eq(2), X3 - current_eq(3),P);
figure(2)
mintgreen = [0.5 1 0.5];
hold on
isosurface(X1, X2, X3, V, level_candidate)
scatter3(ellipseBoxVerts(:,1), ellipseBoxVerts(:,2),ellipseBoxVerts(:,3), 50, mintgreen, 'filled');
xlabel('angle')
ylabel('angular velocity')
zlabel('torque')
hold off
title('Ellipsoid target sets and the vertices of their bounding boxes')
end
%% Compute toe target set
eq_angle = theta_ub(1)-0.03;

current_eq = [eq_angle, 0, -m*g*l*sin(eq_angle) + motor_lim]';

% Linearized system with saturated exo
n = 3;
A = [0, 1, 0; (g/l)*cos(eq_angle), 0, 1/(m*l*l); 0, 0, 0];
B = [0;0;1];


% SDP
Q = sdpvar(n,n); 

Y = sdpvar(1,n);


eps = 1; 

Consts = [[Q>=eps*eye(n)]; [Y'*B'+Q*A' + A*Q + B*Y]<=0*eye(n)];

Obj = ([trace(Q)]); % multiply Q with a scaling diagonal

ops = sdpsettings('solver', 'sedumi');

diagnostics = solvesdp(Consts, Obj, ops);

if diagnostics.problem == 0

    disp('Feasible: CLF found')
    Qval = double(Q); %use "double" to extract solution
    Yval = double(Y);
    P = inv(Qval);
    K = Yval*P

end
R = chol(P);
%% Find level
sway_velocity = 0.01;
% TODO - torque ranges can be found from CoP lb and ub constraint values
% Construct constraint polytope 
A_con = [K; -K; 1 0 0; -1 0 0; 0 1 0;0 -1 0;0 0 1;0 0 -1];
b_con = [uMax; -uMin; theta_ub; -theta_lb; sway_velocity; sway_velocity; tau_max; -tau_min;] + A_con*current_eq;
verts = lcon2vert(A_con,b_con)


% select level set candidate
level_candidate = computeLevelCandidate(A_con, b_con, Qval, R, current_eq);


[axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate)

% convert to vertices
ellipseBoxVerts = boxVerts3d(axes_max, current_eq);
size(ellipseBoxVerts)

% convert vertices to linear constraints using vert2lcon
[A_ts, b_ts] = vert2lcon(ellipseBoxVerts);

% Compute the vertices of the constraint polytope intersected with the
% bounding box of the max inscribed elliposid
verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts]);

target_tol = 0.0001;

% If the intersection does not have 8 vertices, it is not fully contained
% inside the constraint polytope

while size(verticesIntersection,1) ~= 8  

    level_candidate = level_candidate*0.95;
    [axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate);
    ellipseBoxVerts = boxVerts3d(axes_max, current_eq);
    [A_ts, b_ts] = vert2lcon(ellipseBoxVerts);
    verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts]);

end


while max(abs(verticesIntersection - ellipseBoxVerts),[], 'all') > target_tol 
    level_candidate = level_candidate*0.95;
    [axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate);
    ellipseBoxVerts = boxVerts3d(axes_max, current_eq);
    [A_ts, b_ts] = vert2lcon(ellipseBoxVerts);
    verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts]);
end
if plotting
figure(2)
f = @(x1, x2, x3, S) S(1,1).*x1.*x1 + S(1,2).*x1.*x2 + S(1,3).*x1.*x3 + S(2,1).*x1.*x2 + S(2,2).*x2.*x2 + S(2,3).*(x3.*x2)...
    + S(3,3).*x3.*x3  + S(3,1).*x1.*x3 + S(3,2).*x2.*x3;
V = f(X1 - current_eq(1), X2 - current_eq(2), X3 - current_eq(3),P);
mintgreen = [0.5 1 0.5];
hold on
isosurface(X1, X2, X3, V, level_candidate)
scatter3(ellipseBoxVerts(:,1), ellipseBoxVerts(:,2),ellipseBoxVerts(:,3), 50, mintgreen, 'filled');
xlabel('angle')
ylabel('angular velocity')
zlabel('torque')
hold off
title('Ellipsoid target sets and the vertices of their bounding boxes')
end
theta_lb = theta_lb + pi/2;
theta_ub = theta_ub + pi/2;
TS_array = [TS_array (ellipseBoxVerts - current_eq')];

%% Save TS bounds
targetSets.thetaAnkRange = abs(TS_array(1,1));
targetSets.thetaToeRange = abs(TS_array(1,4));
targetSets.thetaDotAnkRange = abs(TS_array(1,2));
targetSets.thetaDotToeRange = abs(TS_array(1,5));
targetSets.torqueAnkRange = abs(TS_array(1,3));
targetSets.torqueToeRange = abs(TS_array(1,6));

%% Pack Target set and corresponding problem parameters
if saveData

    % targetSets, exoParams, alphaMT, alphaRTD
    problemParams.targetSets = targetSets;
    problemParams.thetaLb = theta_lb;
    problemParams.thetaUb = theta_ub;
    
    % exo params
    exoParams.Kp = Kp;
    exoParams.Kd = Kd;
    exoParams.muGrav = muGrav;
    exoParams.motorLim = motor_lim;

    problemParams.exoParams = exoParams;
    
    % Scaling params
    problemParams.alphaMT = alphaMT;
    problemParams.alphaRTD = alphaRTD;
    % Params_YF_noExo_alphaMT_1_alphaRTD_1
    filename = strcat('Params_', modelLabel,'_',exoLabel,'_', num2str(alphaMT*100),'_pctMT_', num2str(alphaRTD*100),'_pctRTD', '.mat')
    save(filename, "problemParams");
end

end