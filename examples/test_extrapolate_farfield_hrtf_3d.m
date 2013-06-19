%TEST_EXTRAPOLATE_FARFIELD_HRTF_3D tests the far field range extrapolation
% of a given spherical HRTF dataset. It is necessary to have the dataset as 
% irs-structure in the folder where this function is located.

%% ===== Configuration ===================================================
clear all; close all; clc;

% load the configuration
conf = SFS_config;

% add to path
addirspath;

% defining grid 
conf.array = 'spherical';
conf.grid = 'HRTFgrid'; 
conf.dimension = '3D';

% load HRIR data set
irs = read_irs('FABIAN_3D_anechoic.mat');
% load('irs_fs.mat');
% irs = irs_fs;

% distance of ps / fs from the listener position
% r = 0.5;

%% ===== Computation =====================================================
% far field extrapolation
far_field_irs = extrapolate_farfield_hrtfset3d(irs,conf);
% near field extrapolation
% near_field_irs = extrapolate_nearfield_hrtfset3d(irs,r,conf);
