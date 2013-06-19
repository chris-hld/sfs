function far_field_irs = extrapolate_farfield_hrtfset_3d(irs,conf)
%EXTRAPOLATE_FARFIELD_HRTFSET_3D far-field extrapolation of a given 3D HRTF 
%                                dataset
%
%   Usage: far_field_irs = extrapolate_farfield_hrtfset(irs,[conf])
%
%   Input parameters:
%       irs     - IR data set for the virtual secondary sources
%                 the irs dataset has to be stored as defined in 
%                 https://dev.qu.tu-berlin.de/projects/measurements/wiki/IRs_file_format
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       far_field_irs  - IR data set extra polated to conation plane wave IRs
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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
% check_irs(irs);
if nargin<nargmax
    conf = SFS_config;
else
%     isargstruct(conf);
end


%% ===== Configuration ===================================================
fs = conf.fs;                   % sampling frequency


%% ===== Variables ======================================================
conf.array = 'spherical';
conf.usetapwin = 0;
conf.usehpre = 0;
conf.usefracdelay = 0;
conf.xref = [0 0 0];

% define secondary source positions and calculate the direction vectors
x0 = zeros(length(irs.source_position),6);
x0(:,1:3) = irs.source_position.';
x0(:,4:6) = direction_vector(x0(:,1:3),repmat(conf.xref,length(irs.source_position),1));

% calculate aliasing frequency for the used grid << approximation >>
conf.hprefhigh = aliasing_frequency_3d(x0(:,1:3));
conf.hpreflow = 1;

%% ===== Computation =====================================================
% get virtual secondary source positions
x0_all = x0;
% Initialize new irs set
far_field_irs = irs;
far_field_irs.description = 'Extrapolated HRTF set containing plane waves';
far_field_irs.left = zeros(size(far_field_irs.left));
far_field_irs.right = zeros(size(far_field_irs.right));
far_field_irs.distance = 'Inf';

% Generate a irs set for all given angles
for ii = 1:length(irs.apparent_azimuth)
    disp([ii length(irs.apparent_azimuth)])

    % direction of plane wave
    xs = -[cos(irs.apparent_azimuth(ii)).*cos(irs.apparent_elevation(ii)) ...
           sin(irs.apparent_azimuth(ii)).*cos(irs.apparent_elevation(ii)) ...
           sin(irs.apparent_elevation(ii))];
    % calculate active virtual speakers
    x0 = secondary_source_selection(x0_all,xs,'pw');

    % calculate driving function
    [~,delay,weight] = driving_function_imp_wfs(x0,xs,'pw',conf);
        
    % sum up contributions from individual virtual speakers
   for l=1:size(x0,1)
        % get IR for the secondary source position
        [phi,theta,r] = cart2sph(x0(l,1),x0(l,2),x0(l,3));
        dt = delay(l)*fs;
        ir_tmp = get_ir(irs,phi,theta,r,[0 0 0]);
        % truncate IR length
        irl = fix_ir_length(ir_tmp(:,1),length(ir_tmp(:,1)),0);
        irr = fix_ir_length(ir_tmp(:,2),length(ir_tmp(:,2)),0);
        % delay and weight HRTFs
        far_field_irs.left(:,ii) = far_field_irs.left(:,ii) + delayline(irl',dt,weight(l),conf)';
        far_field_irs.right(:,ii) = far_field_irs.right(:,ii) + delayline(irr',dt,weight(l),conf)';
    end

end

%% ===== Pre-equalization ===============================================
conf.usehpre = 1;
far_field_irs.left = wfs_preequalization3d(far_field_irs.left,conf);
far_field_irs.right = wfs_preequalization3d(far_field_irs.right,conf);
