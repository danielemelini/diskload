function [u,v,g]=diskload(alpha,icomp,theta,w,nmin,nmax,h,k,l)
%DISKLOAD elastic response to a uniform circular load on a spherical earth
% This function computes the response to a uniform surface pressure load 
% imposed in a disc of angular radius alpha. The elastic response is found 
% at one or more stations located on the surface of the earth at angular 
% distance(s) theta from the center of the disc load.
% The elastic response is computed using user-supplied elastic loading 
% Love numbers (h,k,l) generated using a specific elastic structure model
% for the earth. If three output arguments are invoked, this function
% also computes the change in the height of the geoid at each station.
% The pressure load imposed within the disk is expressed in terms of the 
% equivalent depth height, thickness) of liquid water (density=1000 kg/m3).
%
% USAGE:
%   [u,v] = diskload(alpha,icomp,theta,w,nmin,nmax,h,k,l);
% [u,v,g] = diskload(alpha,icomp,theta,w,nmin,nmax,h,k,l);
%
% INPUT VARIABLES:
%   alpha   disc half-amplitude (degrees)
%   icomp   switch for a compensated (1) / uncompensated (0) disc
%               load
%   theta   angular distances of stations from disc center (degrees)
%       w   pressure imposed on the surface within the spherical disk
%           expressed as the height or depth of equivalent water load (m)
%    nmin   minimum harmonic degree of the expansion to be used
%    nmax   maximum harmonic degrees of the expansion to be used
%             (may be a scalar or a vector with multiple truncation points)
%       h   (n+1)-vector containing the loading Love number h for degrees 0:n
%       k   (n+1)-vector containing the loading Love number k for degrees 0:n
%       l   (n+1)-vector containing the loading Love number l for degrees 0:n
%
% OUTPUT VARIABLES:
%     u     radial or 'vertical' elastic displacement (mm)
%     v     tangential or 'horizontal' elastic displacement (mm)
%     g     geoid change (mm)
%     
%     Size of output vectors is [ length(nmax) length(theta) ]
%     e.g.,   u(i,:) represents U vs theta at the i-th value of nmax
%             u(:,j) represents U vs nmax at the j-th value of theta
%
% NOTES:
% (1) All elements of nmax must be <= n, the maximum order provided for
%     the elastic loading Love numbers
% (2) input w can be positive or negative allowing the user to model
%     the incremental response to incremental pressure changes
% (3) It is easy to switch from the state to the rate problem. If input
%     h is actually the rate of change of the load (in m/yr w.e.), then 
%     outputs u,v and g will become velocities (in mm/yr).
%  
% REFERENCE:
% This function is associated with the publication 
%    Bevis, M., Melini, D., Spada, G., 2016. On computing the geoelastic 
%    response to a disk load, Geophys. J. Int. 205 (3), 1,804-1,812,
%    doi:10.1093/gji/ggw115
% 
% DEPENDENCIES:
% This function calls function pLegendre.m
%
 
% v 1.0 DM 03.07.2015  -- original version ported from REAR
% v 1.1 MB, 
%      reordered LNs in input argument list, to follow the norm
%      switched density of reference load to that of pure water 
%      modified some variable names, changed function name
%      added some input checks
% v 1.2 MB  switch the degree index from l to n so as (i) to aviod
% a conflict with using l as one of the Love numbers, and (ii) follow
% the symbology of the paper. switch the load amplitude parameter from
% h to w, (i) to aviod a conflict with using h as one of the Love numbers,
% and because it reminds us we are using water equivalent representations.
% This variant has function name sdiskload.  s reminds us we are talking of
% a spherical disk.
% v 1.3 DM  replaced the loop over theta with vectorized expressions
% v 1.4 DM  Added degree-0 
% v 1.5 DM  Added support for multiple nmax values

% Define or compute some constants

ggg=6.67384e-11;        % Newton's constant (SI units) 
radius=6371;            % Radius of the Earth (km)
radiusm=radius*1e3;     % Radius of the Earth (m)
grav=9.8046961;         % Surface gravity (m/s/s) 
rhow = 1000;           % density of pure water(kg/m^3) 
rhoear=3.0*grav/4.0/ggg/pi/radiusm;       % Average Earth density (kg/m^3)
from_m_to_mm = 1000;

% check for illegal conditions
m=length(h);
if length(k)~=m || length(l)~=m
    error('Love number vectors do not have same length')
end
if nmax>(m-1)
    error('nmax exceeds the lengths of the Love Number vectors')
end

% check for 0-order LNs
if( (l(1)*k(1)) ~= 0 )
    error('n=0 Love numbers l_0 and k_0 must be zero');
end

% Convert theta to a row vector

dim = size(theta);
if ( dim(2)==1 ), theta = theta'; end

% Computing the harmonic coefficients of the load 
% Vectors are "offset-indexed", i.e.
% P_n = leg(n+1), sigma_n = sigma(n+1)

leg   = pLegendre( max(nmax)+1, cosd(alpha) );
sigma = nan([max(nmax)+1 1]);

switch icomp
    case 0        % Uncompensated disc load, eq. (7)
        for n=nmin:max(nmax)
            if( n==0 ), sigma(n+1) = 0.5*(1-cosd(alpha)); end
            if( n>0  ), sigma(n+1) = 0.5*(-leg(n+2)+leg(n)); end
        end
    case 1        % Compensated disc load, eq. (8)
        for n=nmin:max(nmax)
            if( n==0 ), sigma(n+1) = 0; end
            if( n>0  ), sigma(n+1) = -(leg(n+2)-leg(n)) ... 
                                                / (1+cosd(alpha)); end
        end
end

u = zeros([ length(nmax) length(theta) ]);
v = zeros([ length(nmax) length(theta) ]);
if nargout>2
    g = zeros([ length(nmax) length(theta) ]);
end

% Compute Legendre polynomials at cos(theta)

x   = cosd(theta);
leg = pLegendre( max(nmax)+1, x );

idx = (abs(x) == 1);

% Add the n=0 terms, if required

if( nmin==0 )
    u = u + h(1) * sigma(1);
    if nargout>2
       g = g + sigma(1);
    end
end

for n=max(1,nmin):max(nmax)
       
    dleg= - (n+1) * (x.*leg(n+1,:)-leg(n+2,:)) ...
        ./((1-x).*(1+x)).*sqrt(1-x.^2);
    
    dleg(idx) = 0.;
    
    ampl = sigma(n+1)/(2*n+1);
    
    for i=1:length(nmax)
        if( n<=nmax(i) )
            u(i,:)  =   u(i,:) +     h(n+1)  * ampl * leg(n+1,:);
            v(i,:)  =   v(i,:) +     l(n+1)  * ampl * dleg;
            if nargout>2
                g(i,:)  =   g(i,:) +  (1+k(n+1)) * ampl * leg(n+1,:);
            end
        end
    end    
end

u = u * (3*rhow/rhoear) * w * from_m_to_mm;
v = v * (3*rhow/rhoear) * w * from_m_to_mm;
if nargout>2
    g = g * (3*rhow/rhoear) * w * from_m_to_mm;
end

end

    

        
        

