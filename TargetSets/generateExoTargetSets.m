function [TS_array] = generateExoTargetSets(model,alphaRTD, alphaMT, muGrav, motor_lim)
%"generateExoTargetSets" Function to compute box-shaped target sets
%
%   Function takes in a model that includes mass, height, and bounds on RTD
%   and torque. The user should supply the desired scaling factors alphaRTD and
%   alphaMT and the percentage of gravity compenstion muGrav, and the exoskeleton
%   motor limit.
%
%   Function output is an 8x2 array of box coordinates centered at the
%   ankle and toe equilibrium points.

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
else
    Kp = m*g*l;
    Kd = 0.3*sqrt(m*l*l*Kp);
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

%% Adjust theta_lb and ub if there is enough torque
theta_lb = max(theta_lb, theta_min - (pi/2));
theta_ub = min(theta_ub, theta_max - (pi/2));

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
% chol(P)

%% Find level
sway_velocity = 0.01;

% Construct constraint polytope 
A_con = [K; -K; 1 0 0; -1 0 0; 0 1 0;0 -1 0;0 0 1;0 0 -1];
b_con = [uMax; -uMin; theta_ub; -theta_lb; sway_velocity; sway_velocity; tau_max; -tau_min;]
verts = lcon2vert(A_con,b_con);

% Select level set candidate (max ellipsoid contained in the constraints)
level_candidate = computeLevelCandidate(A_con, b_con, Qval);


% Compute the minimum circumscribed box 
[axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate);

% convert to vertices
ellipseBoxVerts = boxVerts3d(axes_max, current_eq);

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

%% Find level
sway_velocity = 0.01;
% TODO - torque ranges can be found from CoP lb and ub constraint values
% Construct constraint polytope 
A_con = [K; -K; 1 0 0; -1 0 0; 0 1 0;0 -1 0;0 0 1;0 0 -1];
b_con = [uMax; -uMin; theta_ub; -theta_lb; sway_velocity; sway_velocity; tau_max; -tau_min;]
verts = lcon2vert(A_con,b_con)


% select level set candidate
level_candidate = computeLevelCandidate(A_con, b_con, Qval);


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

TS_array = [TS_array (ellipseBoxVerts - current_eq')];

end