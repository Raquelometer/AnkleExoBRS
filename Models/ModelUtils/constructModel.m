function [] = constructModel(label, m, h, tauMin,tauMax, uMin, uMax)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
model.m = m;
model.h = h;
model.tauMin = tauMin;
model.tauMax = tauMax;
model.uMin = uMin;
model.uMax = uMax;
model.label = label;
filename = strcat(label, '_model.mat')
save(filename, "model")
end