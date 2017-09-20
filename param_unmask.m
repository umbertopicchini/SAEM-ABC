
function bigtheta = param_unmask(theta,PARMASK,PARBASE)
%
%         Returns the complete parameter set by inserting the running theta free subset parameter values 
%         into the saved PARBASE all parameter values. 
%
% usage:  bigtheta = param_unmask(theta,parmask,parbase)
% IN:     theta (n,1); running theta free subset parameter
%         parmask; vector 0/1 specifying which parameter are free (1) and which one are not-free (0)
%         parbase; all saved parameter values
% OUT:    bigtheta (m,1); complete vector of parameters updating with the running theta free parameters

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

bigtheta = PARBASE;
% filter theta into bigtheta according to PARMASK:
bigtheta(PARMASK>0) = theta;


