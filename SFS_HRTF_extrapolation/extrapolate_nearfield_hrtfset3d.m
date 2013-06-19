function near_field_irs = extrapolate_nearfield_hrtfset_3d(irs,r,conf)
%EXTRAPOLATE_FARFIELD_HRTFSET_3D far-field extrapolation of a given 3D HRTF 
%                                dataset
%
%   Usage: near_field_irs = extrapolate_nierfield_hrtfset_3d(irs,r,[conf])
%
%   Input parameters:
%       irs     - IR data set for the virtual secondary sources
%                 the irs dataset has to be stored as defined in 
%                 https://dev.qu.tu-berlin.de/projects/measurements/wiki/IRs_file_format
%       conf    - optional configuration struct (see SFS_config)
%       r       - distance of the position of the point source / focused
%                 source to the listener position
% 
%   Output parameters:
%       near_field_irs  - IR data set extra polated 
%
%   EXTRAPOLATE_FARFIELD_HRTFSET_3D(IRS) generates a far-field extrapolated 
%   set of impulse responses, using the given irs set. Far-field means that 
%   the resulting impulse responses are plane waves. The extrapolation is 
%   done via 3D WFS.
%
%   see also: ir_point_source, get_ir, driving_function_imp_wfs_3d

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut für Nachrichtentechnik                    *
%                         Universität Rostock                                *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
end

%% ===== Configuration ===================================================
fs = conf.fs;                   % sampling frequency


%% ===== Variables ======================================================
conf.array = 'spherical';
conf.usetapwin = 0;
conf.usehpre = 0;
conf.usefracdelay = 0;
conf.xref = [0 0 0];

% load irs source positions from HRTF dataset 
x0 = zeros(length(irs.source_position),7);
x0(:,1:3) = irs.source_position.';
% calculate direction vectors for source positions
x0(:,4:6) = direction_vector(x0(:,1:3),repmat(conf.xref,length(irs.source_position),1));
% apply surface weights of the sphere
x0(:,7) = (irs.distance.^2)' .* cos(irs.apparent_elevation)';

%% ===== Computation =====================================================
% get virtual secondary source positions
x0_all = x0;
% Initialize new irs set for ps / fs
near_field_irs = irs;
near_field_irs.description = 'Extrapolated near field HRTF set';
near_field_irs.left = zeros(size(near_field_irs.left));
near_field_irs.right = zeros(size(near_field_irs.right));
near_field_irs.distance = r;

% calculate the secondary source positions for the new data set
[x,y,z] = sph2cart(irs.apparent_azimuth(:),irs.apparent_elevation(:),near_field_irs.distance);
near_field_irs.source_position = [x y z].';

% calculate the aliasing frequency for the used grid << approximation >>
conf.hprefhigh = aliasing_frequency_3d(near_field_irs.source_position(:,1:3));
conf.hpreflow = 1;

% Generate a irs set for all given angles
for ii = 1:length(irs.apparent_azimuth)
    disp([ii length(irs.apparent_azimuth)])

   % direction of plane wave
    xs = near_field_irs.distance.*[cos(irs.apparent_azimuth(ii)).*cos(irs.apparent_elevation(ii)), ...
            sin(irs.apparent_azimuth(ii)).*cos(irs.apparent_elevation(ii)), ...
            sin(irs.apparent_elevation(ii))];

    x0 = secondary_source_selection(x0_all,xs,'fs',conf.xref);

    % calculate driving function
    [~,delay,weight] = driving_function_imp_wfs(x0,xs,'fs',conf);
    
    % sum up contributions from individual virtual speakers
   for l=1:size(x0,1)
        dt = delay(l,1)*fs; 
        % get IR from dataset for the secondary source position
        [phi,theta,r] = cart2sph(x0(l,1),x0(l,2),x0(l,3));
        ir_tmp = get_ir(irs,phi,theta,r,conf.xref');
        % truncate IR length
        ir_tmp = fix_ir_length(ir_tmp,length(ir_tmp(:,1)),dt);
        irl = ir_tmp(:,1);
        irr = ir_tmp(:,2);
        % delay and weight HRTFs
        near_field_irs.left(:,ii) = near_field_irs.left(:,ii) + delayline(irl',dt,weight(l),conf)';
        near_field_irs.right(:,ii) = near_field_irs.right(:,ii) + delayline(irr',dt,weight(l),conf)';
    end

end



%% ===== Pre-equalization ===============================================
conf.usehpre = 1;
near_field_irs.left = wfs_preequalization3d(near_field_irs.left,conf);
near_field_irs.right = wfs_preequalization3d(near_field_irs.right,conf);
