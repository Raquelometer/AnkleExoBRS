function GRFConstraintMat = generateLiftedInvPendExoGRFConstraint(vs, shape, obj)
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


exo_controller = @(x1, x2) max(min((Kd*(0-x2) + Kp*((pi/2)-x1)), motor_lim),-motor_lim) - mu_grav*m*g*l*saturatedCosVal(x1,obj.m,obj.l,motor_lim, 0.01);


c_grf_ub_nom = @(x1, x2, x3) x3 - ((m*l*l)*x2*x2*tan(x1) - ((m+mf)*g*l)/cos(x1) + m*g*l*cos(x1));
c_grf_lb_nom = @(x1, x2, x3) -x3 + ((m*l*l)*x2*x2*tan(x1) - ((m+mf)*g*l)/cos(x1) + m*g*l*cos(x1));

c_grf_ub = @(x1, x2, x3) c_grf_ub_nom(x1,x2,x3) - exo_controller(x1,x2);
c_grf_lb = @(x1,x2,x3) c_grf_lb_nom(x1,x2,x3) + exo_controller(x1,x2);
%GRFConstraintMat = -ones(shape)*Inf;
GRFConstraintMat = ones(shape)*Inf;


for i = 1:shape(1)
    for j = 1:shape(2)
        for k = 1:shape(3)
            if cos(vs{1}(i)) <= 0
                % GRF constraint acts as upper bound
                GRFConstraintMat(i,j,k) = c_grf_lb(vs{1}(i),vs{2}(j),vs{3}(k));
            elseif cos(vs{1}(i)) > 0
                GRFConstraintMat(i,j,k) = c_grf_ub(vs{1}(i),vs{2}(j),vs{3}(k));
            end
        end
    end
end



end