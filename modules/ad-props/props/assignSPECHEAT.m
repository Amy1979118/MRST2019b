function f = assignSPECHEAT(f, specheat, reg)
   % Compute tables (static data)
   kn =  {'uG','uW','uO'}
   for i = 1:3 
     tab = cellfun(@(x)x(:,[1,1+i]), specheat, 'UniformOutput', false);
     %NB the table is of spesific heat the energy is then the integral but we
     % neglect this
     f.(kn{i}) =@(T,varargin) func(T,tab,reg,varargin{:});
   end
end
function v = func(x, tab, reg, varargin)
inx = getRegMap(x, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), tab, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, x, inx);
%NB
v=v.*x;
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
