function [hpreflow, hprefhigh] = localwfs_findhpref(X, phi, xs, src, conf)

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************

%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

if conf.debug
    isargposition(X);
    isargxs(xs);
    isargscalar(phi);
    isargchar(src);
end

%% ===== Configuration ==================================================
fs = conf.fs;               % Sampling rate
conf.N = 2^14;              % Number of Samples
N = conf.N;                 
dimension = conf.dimension; % dimensionality

conf.wfs.usehpre = false;     % no prefilter
conf.localsfs.wfs = conf.wfs;
%% ===== Variables ======================================================
irs = dummy_irs(N);         % Impulse responses

%% ===== Computation ====================================================
% Compute impulse response/amplitude spectrum without prefilter
ir = ir_localwfs(X, phi, xs, src, irs, conf);
[H,~,f]=easyfft(ir(:,1),conf);

H = H./H(1);  % Normalize amplitude spectrum with H(f=0Hz)

if strcmp('2.5D',dimension)
    % Model of local WFS spectrum without prefilter:
    %            _______________
    %           | flow^2+fhigh^2
    %  H(f) = \ |---------------  for f <= fhigh
    %          \|   flow^2+f^2
    %
    %  H(flow),/H(0) = 1/sqrt(2)
    Xfilt = @(x, xlow) sqrt( (xlow.^2)./(xlow.^2 + x.^2));

elseif strcmp('3D',dimension) || strcmp('2D',dimension)
    % Model of local WFS spectrum without prefilter:
    %         
    %          flow^2+fhigh^2
    %  H(f) = ---------------  for f <= fhigh
    %           flow^2+f^2
    %
    %  H(flow)./H(0) = 1/sqrt(2)
    Xfilt = @(x, xlow) (xlow.^2)./(xlow.^2 + x.^2);
else
    error('%s: %s is not a valid conf.dimension entry',upper(mfilename));
end

flowidx = find(H <= Xfilt(1,1), 1, 'first');
hpreflow = f(flowidx);

fhighidx = flowidx - 1 + ...
  find(H(flowidx:end) >= Xfilt(f(flowidx:end), hpreflow), 1, 'first');
  
if isempty(fhighidx)
  hprefhigh = fs/2;
else
  hprefhigh = f(fhighidx);
end 

end
