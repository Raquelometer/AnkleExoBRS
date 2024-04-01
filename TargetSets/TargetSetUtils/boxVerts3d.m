function [verts] = boxVerts3d(axes, eqPt)

% Compute vertices of the minimum circumscribed box
% Shift if eq is not at origin.

vert1 = [-axes(1) -axes(2) -axes(3)] + eqPt';
vert2 = [-axes(1) -axes(2) axes(3)] + eqPt';
vert3 = [-axes(1) axes(2) -axes(3)] + eqPt';
vert4 = [-axes(1) axes(2) axes(3)] + eqPt';
vert5 = [axes(1) -axes(2) -axes(3)] + eqPt';
vert6 = [axes(1) -axes(2) axes(3)] + eqPt';
vert7 = [axes(1) axes(2) -axes(3)] + eqPt';
vert8 = [axes(1) axes(2) axes(3)] + eqPt';


verts =[vert1; vert2; vert3; vert4; vert5; vert6; vert7; vert8;];




end