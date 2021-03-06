function ir = ir_generic(X,phi,x0,d,sofa,conf)
%IR_GENERIC Generate a IR
%
%   Usage: ir = ir_generic(X,phi,x0,d,sofa,conf)
%
%   Input parameters:
%       X       - listener position / m
%       phi     - listener direction [head orientation] / rad
%                 0 means the head is oriented towards the x-axis.
%       x0      - secondary sources [n x 7] / m
%       d       - driving signals [m x n]
%       sofa    - impulse response data set for the secondary sources
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - Impulse response for the desired driving functions (nx2 matrix)
%
%   IR_GENERIC(X,phi,x0,d,sofa,conf) calculates a binaural room impulse
%   response for the given secondary sources and driving signals.
%
%   See also: ir_wfs, ir_nfchoa, ir_point_source, auralize_ir

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 6;
narginchk(nargmin,nargmax);
isargposition(X);
isargscalar(phi);
isargsecondarysource(x0);
isargmatrix(d);
isargstruct(conf);
if size(x0,1)~=size(d,2)
    error(['%s: The number of secondary sources (%i) and driving ', ...
        'signals (%i) does not correspond.'], ...
        upper(mfilename),size(x0,1),size(d,2));
end

%% ===== Configuration ==================================================
N = conf.N; % target length of BRS impulse responses


%% ===== Variables ======================================================
phi = correct_azimuth(phi);


%% ===== BRIR ===========================================================
% Initial values
ir_generic = zeros(N,2);

% Create a BRIR for every single loudspeaker
warning('off','SFS:irs_intpol');
for ii=1:size(x0,1)

    % === Get the desired impulse response.
    % If needed interpolate the given impulse response set and weight, delay the
    % impulse for the correct distance
    ir = get_ir(sofa,X,[phi 0],x0(ii,1:3),'cartesian',conf);

    % === Sum up virtual loudspeakers/HRIRs and add loudspeaker time delay ===
    % Also applying the weights of the secondary sources including integration
    % weights or tapering windows etc.
    ir_generic = ir_generic + fix_length(convolution(ir,d(:,ii)),N).*x0(ii,7);

end
warning('on','SFS:irs_intpol');


%% ===== Headphone compensation =========================================
ir = compensate_headphone(ir_generic,conf);
