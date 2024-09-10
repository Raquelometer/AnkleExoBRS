function CoPConstraintMat = generateLiftedInvPendExoCoPConstraint(vs, shape, obj)
m = obj.m;
l = obj.l;
scaling = true;

g = obj.g;
ank = obj.ank;
lf = obj.lf;
c = obj.c;
b = obj.b;
mf = obj.mf;


%% Scaling

omega = sqrt(l/9.8);
vscale = 1*omega;
mscale = 1;
tauscale = mscale*(omega^2)*(1/(m*l*l));

% Controller params
Kp = obj.Kp;
Kd = obj.Kd;
mu_grav = obj.muGrav;
motor_lim = obj.motorLim;

if scaling

    lf = lf/l;
    ank = ank/l;
    b = b/l;
    c = c/l;
    mf = mf/m;

    m = 1;
    l = 1;
    g = 1;

end

%% Looped version

c_ub = @(x1, x2) ((b*sin(x1)-(ank)*cos(x1))*m*g*cos(x1)-((ank)*sin(x1)+b*cos(x1))*m*l*x2*x2 + (ank)*(mf+m)*g + c*mf*g)./(1 - ((ank)/l)*cos(x1)+(b/l)*sin(x1));
c_lb = @(x1, x2) ((b*sin(x1)+(lf-ank)*cos(x1))*m*g*cos(x1)+((lf-ank)*sin(x1)-b*cos(x1))*m*l*x2*x2 - (lf-ank)*(mf+m)*g + c*mf*g)./(1 + ((lf-ank)/l)*cos(x1)+(b/l)*sin(x1));
exo_controller = @(x1, x2) max(min((Kd*(0-x2) + Kp*((pi/2)-x1)), motor_lim),-motor_lim) - mu_grav*m*g*l*saturatedCosVal(x1,obj.m,obj.l,motor_lim, 0.01);


c_ub_exo = @(x1, x2, x3) -x3 + c_ub(x1,x2) - exo_controller(x1,x2);

c_lb_exo = @(x1, x2, x3) x3 - c_lb(x1, x2) + exo_controller(x1,x2);

%% scale constraints
if scaling
    c_ub_exo = @(x1, x2, x3) -x3 + c_ub(x1,x2) - exo_controller(x1,x2);

    c_lb_exo = @(x1, x2, x3) x3 - c_lb(x1, x2) + exo_controller(x1,x2);
end
constraints = {c_ub_exo, c_lb_exo};


CoPConstraintMat = generateObstacle3(constraints, vs, shape);

%% Vectorized version
% Kx = u_exo_gravComp_PD(xs, obj);
% 
% % TODO double check signs
% CoP_exo_ub = -xs{3} + CoP_constraint_ub(xs, obj) - Kx;
% CoP_exo_lb = xs{3} - CoP_constraint_lb(xs, obj) + Kx;
% 
% % n1 x n2 x n3 matrix
% safety_matrix = CoP_exo_ub - CoP_exo_lb;
% 
% safety_matrix(CoP_ub < uMin) = -700000;
% safety_matrix(CoP_lb > uMax) = -700000;

% Check that range works with torque limits

end