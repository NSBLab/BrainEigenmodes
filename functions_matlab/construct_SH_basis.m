% Construct spherical harmonics basis matrix at specified points.
%
% SYNTAX: [basis] = construct_SH_basis(degree, sphere_points, dl, real_or_complex);
%
% INPUTS:
%
% degree                - maximum degree of spherical harmonics;
% sphere_points         - n sample points on sphere;          [nx3]
% dl                    - {1} for full band; 2 for even order only  
% real_or_complex       - [{'real'} or 'complex'] basis functions.  
%
% OUTPUTS:
%
% basis   - Spherical harmonics values at n points, [ n x k ]
%           where k = (degree+1)^2 for full band or k = nchoosek(degree+2,2) for even only
%
% Example:
%    sphere_points = rand(100,3);
%    [basis] = construct_SH_basis(4, sphere_points, 2, 'real');
function [basis, theta, phi] = construct_SH_basis (degree, sphere_points, dl, real_or_complex)
%%=============================================================
%% Project:   Spherical Harmonics
%% Module:    $RCSfile: construct_SH_basis.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2007/12/27 06:23:35 $
%% Version:   $Revision: 1.8 $
%%=============================================================
if (nargin<2) 
    error('Usage: [basis] = construct_SH_basis(degree(even), points[nx3], [{real}_or_complex]);');  
end
if (nargin<3) 
    dl = 1;  
end
if (nargin<4) 
    real_or_complex = 'real';  
end
n = size(sphere_points,1);  % number of points
[phi, theta] = cart2sph(sphere_points(:,1),sphere_points(:,2),sphere_points(:,3));
%% NOTE: spherical harmonics' convention: PHI - azimuth, THETA - polar angle
theta = pi/2 - theta;  % convert to interval [0, PI]
if dl==1
    k =(degree + 1)*(degree + 1);  %% k is the # of coefficients
else
  if (dl==2)
  	k =(degree + 2)*(degree + 1)/2;  %% k is the # of coefficients
  else
    error('dl can only be 1 or 2');
  end
end
Y = zeros(n,k);
for l = 0:dl:degree   %%% even order only due to the antipodal symmetry!
    % calculate the spherical harmonics
    Pm = legendre(l,cos(theta')); % legendre part
    Pm = Pm';
    lconstant = sqrt((2*l + 1)/(4*pi));
    if (dl==2) 
        center = (l+1)*(l+2)/2 - l;
    else
        if (dl==1) 
            center = (l+1)*(l+1) - l;
        end
    end
    Y(:,center) = lconstant*Pm(:,1);
    for m=1:l
        precoeff = lconstant * sqrt(factorial(l - m)/factorial(l + m));
        switch lower(real_or_complex)
            case 'real'
                 if mod(m,2) == 1
                     precoeff = -precoeff;
                 end
                Y(:, center + m) = sqrt(2)*precoeff*Pm(:,m+1).*cos(m*phi);
                Y(:, center - m) = sqrt(2)*precoeff*Pm(:,m+1).*sin(m*phi);
            case 'complex'
                if mod(m,2) == 1
                    Y(:, center + m) = precoeff*Pm(:,m+1).*exp(i*m*phi);
                    Y(:, center - m) = -precoeff*Pm(:,m+1).*exp(-i*m*phi);
                else
                    Y(:, center + m) = precoeff*Pm(:,m+1).*exp(i*m*phi);
                    Y(:, center - m) = precoeff*Pm(:,m+1).*exp(-i*m*phi);
                end
            otherwise
                error('The last argument must be either \"real\" (default) or \"complex\".');
        end
    end
end
basis = Y;