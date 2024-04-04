function [TS_array, theta_lb, theta_ub] = generateNominalTargetSets(model,modelLabel, alphaRTD, alphaMT)

%"generateNominalTargetSets" Function to compute box-shaped target sets
%
%   Function takes in a model that includes mass, height, and bounds on RTD
%   and torque. The user should supply the desired scaling factors alphaRTD and
%   alphaMT.
%
%   Function output is an 8x2 array of box coordinates centered at the ankle
%   and toe eq points, and the lower and upper boundary of feasible static support 
%   (i.e. over which range of angles can the model maintain a statically stable configuration)
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


%% Compute boundary of foot

theta_heel = acos(ank/l);
theta_toe = pi - acos((lf-ank)/l);

heel_torque = m*9.8*l*cos(theta_heel);

toe_torque = m*9.8*l*cos(theta_toe);

%% Compute static support boundaries

theta_lb = acos(tau_max/(m*g*l)) - (pi/2);
theta_ub = acos(tau_min/(m*g*l)) - (pi/2);

%% Plot relevant values as a sanity check

theta_vals = linspace(0, pi, 100000);
Tau_vals = m*g*l*cos(theta_vals);

figure(1)
hold on
plot(theta_vals - (pi/2), Tau_vals);
yline(tau_min, 'Label', 'Torque lower bound')
yline(tau_max, 'Label', 'Torque upper bound')
xline(theta_heel - (pi/2), 'Color','red','LineWidth',1.5, 'Label','Heel')
xline(theta_toe - (pi/2), 'Color','blue','LineWidth',1.5,'Label','Toe')
xline(theta_lb,'Color','green','LineWidth',1.5, 'Label','Static support lower bound', 'LabelVerticalAlignment','middle')
xline(theta_ub, 'Color','green','LineWidth',1.5, 'Label','Static support upper bound', 'LabelVerticalAlignment','middle')
%xline(pi/2, 'Color', 'green','LineWidth',1)
xlim([-0.2 .5]);
ylim([-150 50]);
hold off

%% Step 2: Adjust theta_lb and theta_ub in case there is more than enough torque
theta_lb = max(theta_lb, theta_heel - (pi/2));
theta_ub = min(theta_ub, theta_toe - (pi/2));

%% Pack eq points
eq_angles = [0; theta_ub - 0.03;] % ankle and toe equilibria

TS_array = zeros(8, 3*length(eq_angles));

%% Set up 3d Plots
if plotting

    x1 = linspace(-1,1, 300);
    x2 = linspace(-1,1,300);
    x3 = linspace(-150,50,500);
    
    [X1, X2, X3] = meshgrid(x1, x2, x3);
end
%% Compute Target sets for each eq point
for i = 1:length(eq_angles)
eq_angle = eq_angles(i);
targetSets.theta_toe_eq = eq_angle + pi/2;
current_eq = [eq_angle, 0, -m*g*l*sin(eq_angle)]';

% Linearized system without exo
n = 3;
A = [0, 1, 0; (g/l)*cos(eq_angle), 0, 1/(m*l*l); 0, 0, 0];
% A = [0, 1, 0; (g/l)*cos(0), 0, 1/(m*l*l); 0, 0, 0];

B = [0;0;1];

% SDP
Q = sdpvar(n,n); 

Y = sdpvar(1,n);

eps = 1; 

Consts = [[Q>=eps*eye(n)]; [Y'*B'+Q*A' + A*Q + B*Y]<=0*eye(n)];

Obj = ([trace(Q)]); 
% Obj = ([]); 


ops = sdpsettings('solver', 'sedumi');

diagnostics = solvesdp(Consts, Obj, ops);

if diagnostics.problem == 0

    disp('Feasible: CLF found')
    Qval = double(Q); %use "double" to extract solution
    Yval = double(Y);
    P = inv(Qval);
    K = Yval*P;
    
end
targetSets.K_heel = K;

R = chol(P);
%% Construct constraint polytope

sway_velocity = 0.01;



A_con = [K; -K; 1 0 0; -1 0 0; 0 1 0;0 -1 0;0 0 1;0 0 -1];
% b_con = [uMax - K_offset; -uMin + K_offset; theta_ub; -theta_lb; sway_velocity; sway_velocity; tau_max; -tau_min;] - A_con*current_eq;
b_con = [uMax; -uMin; theta_ub; -theta_lb; sway_velocity; sway_velocity; tau_max; -tau_min;] + A_con*current_eq;

verts = lcon2vert(A_con,b_con);

%% Compute max inscribed ellipsoid

level_candidate = computeLevelCandidate(A_con, b_con, Qval, R, current_eq)


% Enclose ellipsoid {x:x'Px < level_candidate} inside a minimum circumscribed box
[axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate)

% convert box axes to vertices
ellipseBoxVerts = boxVerts3d(axes_max, current_eq)

% convert target set bounding box vertices back to linear constraints using vert2lcon
[A_ts, b_ts] = vert2lcon(ellipseBoxVerts);

% Compute the vertices of the constraint polytope intersected with the
% bounding box of the max inscribed elliposid
verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts])

target_tol = 0.0001;

level_scaling = 1;

% If the intersection does not have 8 vertices, it is not fully contained
% inside the constraint polytope

while size(verticesIntersection,1) ~= 8  

    level_candidate = level_candidate*0.95;
    [axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate);
    ellipseBoxVerts = boxVerts3d(axes_max, current_eq);
    [A_ts, b_ts] = vert2lcon(ellipseBoxVerts, current_eq);
    verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts]);

    level_scaling = level_scaling + 1;


end


% Shrink TS until the intersection of the constraints and the bounding box are the
% same. NOTE there is an assumption here that the polytopes toolbox uses a
% consistent ordering, which we try to match in the function
% MinCircBox_ellipsoid().

while max(abs(verticesIntersection - ellipseBoxVerts),[], 'all') > target_tol

    level_candidate = level_candidate*0.95;
    [axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level_candidate);
    ellipseBoxVerts = boxVerts3d(axes_max, current_eq);
    [A_ts, b_ts] = vert2lcon(ellipseBoxVerts);
    verticesIntersection = lcon2vert([A_con; A_ts],[b_con; b_ts]);

end

% TS_array = ellipseBoxVerts;

TS_array(:,(i + 2*(i-1)):(i + 2*i)) = ellipseBoxVerts - current_eq';

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

end

theta_lb = theta_lb + pi/2;
theta_ub = theta_ub + pi/2;

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
    exoParams.Kp = 0;
    exoParams.Kd = 0;
    exoParams.muGrav = 0;
    exoParams.motorLim = 0;

    problemParams.exoParams = exoParams;
    
    % Scaling params
    problemParams.alphaMT = alphaMT;
    problemParams.alphRTD = alphaRTD;
    % Params_YF_noExo_alphaMT_1_alphaRTD_1
    filename = strcat('Params_', modelLabel,'_','noExo_', num2str(alphaMT*100),'_pctMT_', num2str(alphaRTD*100),'_pctRTD', '.mat')
    save(filename, "problemParams");
end
end