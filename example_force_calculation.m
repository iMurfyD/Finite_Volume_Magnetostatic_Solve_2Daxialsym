clear all; close all; clc;
addpath('./solver/');
 
%% Set up parameters
sep = 2.2;   % Center to center grain separation distance, normalized to a
             % There is a singularity as sep -> 2 (perfect contact requires
             %                                     infinitely small grid)
H0 = 477.0;  % Applied magnetic field, A/m
susc = 0.96; % Magnetic susceptibility
a = 1.4e-6;  % Grain radius, meters
n = 100;    % Grid size

%% Sets up a two grain finite volume simulation, and then uses a
%  cartesian -> polar interpolation scheme to integrate the Maxwell
%  stress tensor (neglecting electrical contributions) around one of the
%  grains. 
f = calc_f_two_grain(sep, H0, susc,a, n, n)
