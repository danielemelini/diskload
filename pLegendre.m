function p=pLegendre(lmax, z)
%
%   p = pLegendre(lmax, z)
%
%	This subroutine evalutates all of the unnormalized Legendre polynomials 
%	up to degree lmax. 
%
%	Calling Parameters:
%	
%   OUT:
%       p: A vector of all unnormalized Legendgre polynomials evaluated at 
%          z up to lmax. The lenght is equal to (lmax+1).
%   IN
%       lmax: Maximum degree to compute.
%       z:    Value within [-1, 1], cos(colatitude) or sin(latitude).
%
%	Notes:
%	
%	1.	The integral of Pl**2 over (-1,1) is 2/(2l+1).
%	2.	Values are calculated accoring to the following recursion scheme:
%       P_0(z) = 1.0, P_1(z) = z, and 
%       P_l(z) = (2l-1) * z * P_{l-1}(z) / l - (l-1) * P_{l-2}(z) / l
%
%	Dependencies:	None
%
%	Written by Mark Wieczorek June 2004
%       Modified by Giorgio Spada, 2007
%       Ported to MATLAB by Daniele Melini, 2015
%
%	Original code is Copyright (c) 2005, Mark A. Wieczorek
%	All rights reserved.
%

nz   = length(z);

% Convert z to a row-vector

junk = size(z);

if( ( junk(1) ~= 1 ) && ( junk(2) ~= 1 ) )
    error( 'pLegendre: z must be a scalar or a vector.' );
end

if( junk(2) == 1 )
    z = z';
end

if( lmax<0 )
    error( 'pLegendre: lmax must be greater than or equal to 0.' );
end

if( sum( abs(z)>1 ) > 0 ) 
    error( 'pLegendre: abs(z) must be less than or equal to 1.' );
end
      
paux = nan([lmax+2 nz]);
p    = nan([lmax+1 nz]);

pm2       = 1.;
paux(1,:) = 1.;
      	
pm1       = z;
paux(2,:) = pm1;
      	
for l = 2:lmax
    pl = ((2*l-1) * z .* pm1 - (l-1) * pm2)/l;
    paux(l+1,:) = pl;
    pm2  = pm1;
    pm1  = pl;
end

for j=0:lmax
    p(j+1,:)=paux(j+1,:);
end

end
