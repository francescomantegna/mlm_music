%{
    file_name : test_elec.m
    author : Francesco Mantegna
    institution : NYU
    project : Music&Poetry
    date : 12/01/2019
%}

%% add fieldtrip path

addpath(genpath('/Users/francesco/Documents/MATLAB/fieldtrip-20170414'))
ft_defaults

%% load elec

% load('standard_1020_elec.mat')
load('fran_59CH_elec.mat')

x = elec.pnt(:,1);
y = elec.pnt(:,2);
z = elec.pnt(:,3);

%% generate mesh

xv = linspace(min(x), max(x), 1000);
yv = linspace(min(y), max(y), 1000);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);

%% visualize data

figure('Renderer', 'painters', 'Position', [5 5 1000 400]);clf

subplot(1,2,1)

scatter3(x,y,z)

box off
grid on
set(gca, 'TickDir', 'out')
% rotate3d on
view([-29 22])

% figure;

subplot(1,2,2)

surf(X, Y, Z);

box off
grid on
set(gca, 'TickDir', 'out')
% rotate3d on
view([-29 22])

set(gca, 'ZLim',[0 100])
shading interp
