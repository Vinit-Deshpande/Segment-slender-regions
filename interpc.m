function S = interpc(Q,Ns,METHOD)
% Interpolation of a 2D or 3D curve q with simple parameterisation. 
% By John S. Kearns 
% Verions 1.0.0
%__________________________________________________________________________
% Inputs: q - an [Nqx2] or [Nqx3] curve that you wish to interpolate.
%             where Nq = length(q);
%
%         Ns - Number of points desired for interpolation 
%              Default is Ns = 100;
%
%         METHOD - Desired interpolation method, input as char or string.
%                  Default is METHOD = 'makima';
%__________________________________________________________________________
% Output: s - an [(Ns)x2] or [(Ns)x3] curve that is interpolated from q.
%__________________________________________________________________________

%__________________________________________________________________________
% Interpolated curve 'S' will have more dense point spacing in 
% high curvature regions of 'Q'. 'Q' with uniform curvature will produce
% uniform point spacing with arclength in 'S'.
%
% Any method that is available in MATLAB's inbuilt 'interp1()' function is 
% available with 'interpc()'. The default interpolation method is 'makima.'
%
% Code is pretty simple. It is intended to help beginners who might be 
% stumbling their way through MATLAB like I once was. 
% 
% For densely packed input curves 'Q', all interpolation methods tend to 
% give similar results (see test file). For sparse input data, methods 
% differ appreciably in regions of high curvature. See MATLAB documentation 
% on 'interp1()' for details about the behaviour of various interpolation
% methods. 
% 
% Documentation is accessed by entering 
% 
% doc interp1
% 
% in command window.
%
% By John S. Kearns
% Email: john.kearns@monash.edu
%
%   Version #   | Date(Y/M/D) |  Release Notes
%-----------------------------------------------
% Version 1.0.0 |  23/07/21   | Initial Release.
%__________________________________________________________________________
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Notes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function is pretty straightforward. I just wanted to make something
% simple that could easily and robustly interpolate an arbitrary curve.
% Planned Updates: 
% 1) Additional input argument 'even' will assert that output
%    curve s has points of equal arclength spacing. 
%
% 2) Additional interplation options that are not in MATLAB's 'interp1()'.
%
% 3) Accept mixed inputs for METHOD for x, y, & z curve components.  
%
% Function inspired by 'curvspace()', written by Yo Fukushima
% Curvspace Link: 
% https://au.mathworks.com/matlabcentral/fileexchange/...
%                                            7233-curvspace?s_tid=srchtitle
%
% I am jealous of Yo for taking the best name for curve interpolation.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if nargin<3
    METHOD = 'makima';
    if nargin<2
        Ns = 100;
    end
end

% Parameterize source curve Q and target curve S.
s = linspace(0,1,Ns)';
q = linspace(0,1,size(Q,1))';

% Interpolate components of source curve Q 
% Assumes Q(:,1) = qx;
%         Q(:,2) = qy;
%         Q(:,3) = qz;


sx = interp1(q,Q(:,1),s,METHOD);
sy = interp1(q,Q(:,2),s,METHOD);

% Check if curve is 2D or 3D
if size(Q,2) == 3
    sz = interp1(q,Q(:,3),s,METHOD);
end

% Store as single output curve. 
S = [sx,sy];

% Fill if curve is 3D.
if size(Q,2) == 3
    S(:,3) = sz;
end

end


 