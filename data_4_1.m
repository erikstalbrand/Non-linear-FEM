%% Exercise 4.1 a) data.m

clear;
clc;

load('geometry_E41.mat');

% material parameters
E = 210E3;
v = 0.3;
mpara = [E v];

% geometry
w = 0.6;        % width
h = 0.3;        % height
r = 0.05;       % radius
t = 1;          % thickness

ndof = max(max(edof));
nelm = length(edof);
for i = 1:nelm
ec = [ex(i, :); ey(i, :)];       % element coordinates
end

save('data.mat', 'ec', 'ndof', 'nelm', 'mpara', 'w', 'h', 'r', 't', 'bc', 'edof', 'ex', 'ey');

