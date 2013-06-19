%TEST_WAVEFIELD_3D_WFS is a test file. No input parameters are necessary.
% Within the file different properties can be adjusted, e.g. computation of
% a plane wave ('pw'), point source ('ps') or focused source ('fs'). The 
% listener position can be adjusted by conf.xref. The used grid for 3d WFS
% is a spherical one defined by the measurement points of the HRIRs 
% (measured by FABIAN). The frequency of the wavefield in the frequency 
% domain can be adjusted to arbitrary values by editing f.

%% ===== Configuration ===================================================
clc
clear all
close all
% properties of SFS_conf
conf = SFS_config_example;
conf.zreferenceaxis = 'y'; % set a plot reference axis
conf.useplot = 1; % plot the results = 1 ; otherwise = 0

conf.array = 'spherical'; % array type
conf.usetapwin = 0; % do not use tapering window, it's not needed in the 3D case
conf.plot.loudspeakers = 0; % do not plot loudspeakers in the 3D case, because it's a mess ;)
conf.grid = 'HRTFgrid'; % two grids available: MinimumEnergyPoints and HRTFgrid (if you type anything else then 'MinimumEnergyPoints')
conf.number_of_points_on_sphere = 81^2; % number of points on the sphere, if spherical array is choosen
conf.debug = 1; % debug=1 allows to plot results of different evaluation steps
conf.frame = 0; % set a frame to show the wavefield in the time domain
conf.xref = [0 0 0]; % ps/fs: 'listener position' ; pw:  place where the wavefield is scaled to one

% properties of desired wavefield
xs = [0 2 0]; % position of point source or focus source / inicidence angle of plane wave
r = 1.5; % radius of the sphere
L = 2.*r; % diameter of the sphere

x0 = secondary_source_positions(L,conf);

conf.usehpre = 1; % use preequalization filter
conf.hprefhigh = aliasing_frequency_3d(x0(:,1:3),conf); % calculate aliasing frequency << approximation >>
conf.dimension = '3D'; % choose dimension: '2D','2.5D' or '3D'
src = 'ps'; % select source type pw/ps/fs
f = 1000; % frequency at which the wavefield will be calculated
t = 250;  % frame at which the wavefield will be calculated

% plot properties
scale_axis_1 = [0 0];
scale_axis_2 = [-2 2];
scale_axis_3 = [0.5 0.5];
conf.xysamples = 300;

r_str = num2str(r);
point_str_x = num2str(xs(1));
point_str_y = num2str(xs(2));
point_str_z = num2str(xs(3));
f_str = num2str(f);

%% ===== Computation ====================================================
% calculate frequency and time domain wavefields
[x,y,z,ps_weights_f_1000,x0,win] = wave_field_mono_wfs_3d(scale_axis_2,scale_axis_2,scale_axis_1,xs,src,f,L,conf);
title(['Mono-frequent wavefield using 3D WFS with a spherical array, parameters: r = '...
       r_str ', virtual source: ' src ' from/at [' point_str_x ',' point_str_y ',' point_str_z...
       '] and f = ' f_str 'Hz'])

% wave_field_mono_wfs_3d(scale_axis_2,scale_axis_3,scale_axis_2,xs,src,f,L,conf);
% title(['Mono-frequent wavefield using 3D WFS with a spherical array, parameters: r = '...
%        r_str ', virtual source: ' src ' from/at [' point_str_x ',' point_str_y ',' point_str_z...
%        '] and f = ' f_str 'Hz'])

% wave_field_mono_wfs_3d(scale_axis_1,scale_axis_2,scale_axis_2,xs,src,f,L,conf);
% title(['Mono-frequent wavefield using 3D WFS with a spherical array, parameters: r = '...
%        r_str ', virtual source: ' src ' from/at [' point_str_x ',' point_str_y ',' point_str_z...
%        '] and f = ' f_str 'Hz'])

%% calculate wavefield in time domain                                         
% wave_field_imp_wfs(scale_axis_2,scale_axis_2,scale_axis_1,xs,src,t,L,conf);
% grid on
% title('WFS 3D spherical array r = 1.5m, plane wave [0,1,0],time domain')
% wave_field_imp_wfs_3d(scale_axis_2,scale_axis_3,scale_axis_2,xs,src,L,conf);
% grid on
% wave_field_imp_wfs_3d(scale_axis_1,scale_axis_2,scale_axis_2,xs,src,L,conf);
% grid on
%%

