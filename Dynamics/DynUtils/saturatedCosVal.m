function satCos = saturatedCosVal(theta, m, l, motorLim, tWidth)

% Function for computing a smoothed version of the gravity compensation exo
% Not actually necessary.

g = 9.8;

thetaSat_lb = acos(motorLim/(m*g*l)); % theta1
thetaSat_ub = acos(-motorLim/(m*g*l)); % theta2

cosTheta1 = motorLim/(m*g*l);
cosTheta2 = -motorLim/(m*g*l);

satCos = cosTheta1 + 0.5.*(cos(theta) - cosTheta1).*(1 - tanh((theta - thetaSat_lb)/tWidth))...
    + 0.5.*(cosTheta2 - cos(theta)).*(1 - tanh((theta - thetaSat_ub)/tWidth));

end