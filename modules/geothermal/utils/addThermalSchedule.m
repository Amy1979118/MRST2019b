function [ schedule ] = addThermalSchedule( schedule, varargin )
% add thermal boundary condition to an existing schedule structure. 
% 
% SYNOPSIS:
%   schedule = addThermalSchedule(schedule, 'pn1',pv1 ...)
% 
% PARAMETERS:
%   schedule - schedule structure containing info regarding bc, src terms
%              and time steps. Created with simpleSchedule.
%   
%   bcT      - structure that contains the thermal boundary conditions. 
% 
%   srcT     - structure that contains the thermal source terms.
% 
% RETURNS:
%   schedule - valid schedule structure with thermal bc and src terms.
%
% SEE ALSO:
%   'simpleSchedule' , 'makeThermalSchedule'.
 
opt = struct('bcT'  , [] , ... 
             'srcT' , [] );
         
opt = merge_options(opt, varargin{:});

schedule.control.bcT    = opt.bcT; 
schedule.control.srcT   = opt.srcT; 


end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}