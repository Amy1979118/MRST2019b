function [ bc ] = setPressureBC( bc, G, side, type, value, varargin )
% Impose constant pressure or mass flux boundary conditions on global side
% This function combines functions pside and fluxside
% 
% SYNOPSIS 
% bc = setPressureBC(bc, G, side, type, value)
% bc = setPressureBC(bc, G, side, type, value, 'pn', pv)
% bc = setPressureBC(bc, G, side, type, value, I1, I2)
% bc = setPressureBC(bc, G, side, type, value, I1, I2, 'pn', pv)
%
% PARAMETERS
% bc - Boundary condition structure as defined by function 'addBC'
% 
% G  - Grid structure as described by 'grid_structure'. Currently
%      restricted to grids produced with functions 'cartGrid' and 
%      'tensorGrid' and other grids that add cardinal directions to 
%      'G.cells.faces(:,2) in the same format'. 
% 
% side - Global side from which to extract face indices. String. Must 
%       (case insensitively) match one of six alias groups :
%
%               1) `{'West' , 'XMin', 'Left'  }`
%               2) `{'East' , 'XMax', 'Right' }`
%               3) `{'South', 'YMin', 'Back'  }`
%               4) `{'North', 'YMax', 'Front' }`
%               5) `{'Upper', 'ZMin', 'Top'   }`
%               6) `{'Lower', 'ZMax', 'Bottom'}`
%
%            These groups correspond to the cardinal directions mentioned
%            as the first alternative in each group.
%
%   type     - Type of boundary: Neumann - Dirichlet
%              'pressure'   : pressure value in units of Pascal (scalar or vector)    
%              'flux'       : total flux in units of m^3/s (scalar)
%              accounted for by faces on side in ranges I1 and I2 
%
%   value    - Temperature value or flux, to be applied to the face.
%            Either a scalar or a vector of `numel(I1)*numel(I2)` values.
%
%   I1,I2  - Cell index ranges for local (in-plane) axes one and two,
%            respectively.  An empty index range ([]) is interpreted as
%            covering the entire corresponding local axis of 'side' in the
%            grid 'G'.  The local axes on a 'side' in 'G' are ordered
%            according to 'X' before 'Y', and 'Y' before 'Z'.
%
% OPTIONAL PARAMETERS:
%   sat    - Fluid composition of fluid injected across inflow faces.
%            An n-by-m array of fluid compositions with 'n' being the
%            number of individual faces specified by (I1,I2) (i.e.,
%            n==NUMEL(I1)*NUMEL(I2)) or one.  If m=3, the columns of 'sat'
%            are interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is for the benefit of transport solvers such as
%            'blackoilUpwFE' and will be ignored for outflow faces.
%
%            Default value: `sat = []` (assume single-phase flow).
%
%   range  - Restricts the search for outer faces to a subset of the cells
%            in the direction perpendicular to that of the face. Example:
%            if side='LEFT', one will only search for outer faces in the
%            cells with logical indexes [range,:,:].
%            Default value: `range = []` (do not restrict search).
%
% RETURNS:
%   bc     - Updated boundary condition structure.
% 
% SEE ALSO:
% 'pside', 'fluxside', 'setThermalBC' 


% Sanity check
if ~isfield(G,'cartDims'), error(msgid('NotImplemented'),'PressureBC not implemented for this grid type'); end 

mrstNargInCheck(5,10,nargin);

if nargin == 5 || ischar(varargin{1})
    I1 = []; I2 = []; 
else 
    I1 = varargin{1}; I2 = varargin{2};
    varargin = varargin{3:end};
end 

opt = struct('sat',[],'range',[]);
opt = merge_options(opt, varargin{:});
sat = opt.sat;

ix = ThermalboundaryFaceIndices(G, side, I1, I2, opt.range);
    
assert(any(numel(value) == [1, numel(ix)]));
assert (isempty(sat) || any(size(sat,1) == [1, numel(ix)]));


switch type 
    case 'pressure'
        if size(sat,1)  == 1, sat   = sat(ones([numel(ix), 1]), :); end
        if numel(value) == 1, value = value(ones([numel(ix), 1]));  end
        bc = addBC(bc, ix, 'pressure', value, 'sat', sat);
        nf = size(bc.face); 
        bc.hW = zeros(nf);
    case 'flux'
        if size(sat,1)  == 1, sat   = sat(ones([numel(ix), 1]), :); end
        a  = G.faces.areas(ix);
        sa = sum(a);
        bc = addBC(bc, ix, 'flux', (value / sa) .* a, 'sat', sat);
      
    otherwise 
        warning('Unexpected BC type. No BC created');
end 

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

