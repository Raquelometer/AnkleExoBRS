function [axes_min, axes_max] = MinCircBox_ellipsoid(Qval, level)
% MINCIRCBOX_ELLIPSOID computes the tightest axis aligned box about a
% given ellipsoid x'Qval^-1x < level
% Derivation courtesy of Tavian Barnes available at
% https://tavianator.com/2014/ellipsoid_bounding_boxes.html

R = [Qval zeros(3,1); 0 0 0 -1];
% R = [inv(P) zeros(3,1); 0 0 0 -1]
% R(4,4) = -level_ankle
R(4,4) = -level;
zmin = (R(3,4) - sqrt((R(3,4)^2) - (R(4,4)*R(3,3))) ) / R(4,4);
ymin = (R(2,4) - sqrt((R(2,4)^2) - (R(4,4)*R(2,2))) ) / R(4,4);
xmin = (R(1,4) - sqrt((R(1,4)^2) - (R(4,4)*R(1,1))) ) / R(4,4);

zmax = (R(3,4) + sqrt((R(3,4)^2) - (R(4,4)*R(3,3)))) / R(4,4);
ymax = (R(2,4) + sqrt((R(2,4)^2) - (R(4,4)*R(2,2)))) / R(4,4);
xmax = (R(1,4) + sqrt((R(1,4)^2) - (R(4,4)*R(1,1)))) / R(4,4);

% zmax_ankle = zmin*level_ankle
% ymax_ankle = ymin*level_ankle
% xmax_ankle = xmin*level_ankle

zmax = zmin*level;
ymax = ymin*level;
xmax = xmin*level;

axes_max = [xmax; ymax; zmax; ];

axes_min = [-xmax; -ymax; -zmax; ];
end