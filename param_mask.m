function theta = param_mask(bigtheta,PARMASK)

%          Returns the vector of free parameters from the complete vector of
%          parameters bigtheta
%
% usage :  theta = param_mask(bigtheta,parmask)
% IN:      bigtheta (m,1); vector of all parameters
%          parmask; vector 0/1 specifying which parameters are free (1) and which ones are not free (0)
% OUT:     theta (n,1); vector of free parameters

% This file is part of the "abc-sde" program. Copyright (C) 2013 Umberto Picchini
% https://sourceforge.net/projects/abc-sde/
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.



npar = length(PARMASK);

if (npar ~= length(bigtheta))
   error('incompatible dimensions of PARMASK and bigtheta');
end


% interestingly in our case using "find" seems faster than using logical indexing
% theta = bigtheta(PARMASK>0);
theta = bigtheta(find(PARMASK));


