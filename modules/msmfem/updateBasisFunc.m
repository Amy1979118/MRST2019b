function [CS, varargout] = updateBasisFunc(S, CS, G, CG, rock, mob, varargin)
% Update basis functions in regions where the total mobility has changed
% more than a given tolerance.
%
% SYNOPSIS:
%   CS = updateBasisFunc(state, rock, S, CS, CG, mob)
%   CS = updateBasisFunc(state, rock, S, CS, CG, mob'pn1', pv1, ...)
%
% DESCRIPTION:
% Two different modes:
% 1) The user can supply a set of faces where the basis functions should be
% updated.
%
% 2) Update the basis function if any fine cell related to the basis
% function is not within
%
%   1/(1+tol) < mob/mobOld < 1 + tol
%
% where mob is the total mobility and mobOld is the mobility before the
% last saturation update.
%
%
% PARAMETERS:
%   state   - Solution structure.
%
%   rock    - Rock data structure with valid field 'perm'. If the basis
%             functions are to be weighted by porosity, rock must also
%             contain a valid field 'poro'.
%
%   S       - System struture describing the underlying fine grid model,
%             particularly the individual cell flux inner products.
%
%   CS      - Coarse system structure.
%
%   G       - Grid structure as described by grid_structure.
%
%   CG      - Coarse grid structure as defined by generateCoarseGrid.
%
%   mob     - Total mobility.  One scalar value for each cell in the
%             underlying (fine) model.
%
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%               - faces --
%                        Which basis functions to update. Default value:
%                        [], meaning that we use mobOld and mobTol to
%                        calculate what faces to update.
%
%               - mobOld  --
%                        Total mobility from previous time step.  One
%                        scalar value for each cell in the underlying
%                        (fine) model. See mode 2 in function description.
%
%               - mobTol --
%                        Mobility tolerance. See mode 2 in function
%                        description.
%
%               - Verbose --
%                        Whether or not to emit progress reports while
%                        computing basis functions.
%                        Logical.  Default value dependent upon global
%                        verbose setting in function 'mrstVerbose'.
%
%               - bc  -- Boundary condtion structure as defined by function
%                        'addBC'.  This structure accounts for all external
%                        boundary contributions to the reservoir flow.
%                        Default value: bc = [] meaning all external
%                        no-flow (homogeneous Neumann) conditions.
%
%               - src -- Explicit source contributions as defined by
%                        function 'addSource'.
%                        Default value: src = [] meaning no explicit
%                        sources exist in the model.
%
%               - Overlap --
%                        Number of fine-grid cells in each physical
%                        direction with which to extend the supporting
%                        domain of any given basis functions.
%
%                        Using overlapping domains enables capturing more
%                        complex flow patterns, particularly for very
%                        coarse grids, at the expense of increased coupling
%                        in the resulting systems of linear equations.
%                        Non-negative integers.  Default value = 0.
%
%
% RETURNS:
%   CS - System structure with updated fields:
%          - basis   - Flux basis functions as generated by function
%                      evalBasisFunc.
%          - basisP  - Pressure basis functions as generated by function
%
%   faces - List of updated faces.
%
% SEE ALSO:
%   `computeMimeticIP`, `generateCoarseGrid`, `evalBasisFunc`, `mrstVerbose`.

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


[verbose, faces, weight, overlap, src, bc, mobTol, mobOld] = ...
                                    parse_args(G, CS, rock, varargin{:});

if isempty(faces)
   assert(mobTol>=0 & ~isempty(mobOld))
   faces = find_update_faces(CG, CS, mob, mobOld, mobTol);
end


%% Update basis functions
%
cellNo  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
C       = sparse(1:numel(cellNo), cellNo, 1);
D       = sparse(1:numel(cellNo), double(G.cells.faces(:,1)), 1, ...
                 numel(cellNo), G.faces.num);

[V, P] = evalBasisFunc(faces, G, CG, S.BI, C, D,    ...
                       weight, mob, 'src', src, 'bc', bc, ...
                       'Verbose', verbose, 'Overlap', overlap);

CS = assignBasisFuncs(CS, V, P);


if nargout == 2
   varargout{1} = faces;
end

end
%-----------------------------------------------------------------------
% Private helpers follow
%-----------------------------------------------------------------------


function [verbose, faces, weight, overlap, src, bc, mobTol, mobOld] ...
                                       = parse_args(G, CS, rock, varargin)

opt = struct('Verbose',        false,  ...
             'BasisWeighting', 'perm', ...
             'Overlap',        0,      ...
             'src', [], 'bc', [],      ...
             'faces',  [],             ...
             'mobTol', 10,             ...
             'mobOld', [] );

opt = merge_options(opt, varargin{:});

verbose    = opt.Verbose;
weighting  = CS.basisWeighting;
overlap    = opt.Overlap;
src        = opt.src;
bc         = opt.bc;
faces      = opt.faces;
mobTol     = opt.mobTol;
mobOld     = opt.mobOld;

weight     = evalBasisSource(G, weighting, rock);
end

%-----------------------------------------------------------------------

function faces = find_update_faces(CG, CS, mob, mobOld, mobTol)
   cells = (mob./mobOld) > 1+mobTol | (mob./mobOld) <1/(1+mobTol) ;

   coarseBlocks = unique(CG.partition(cells));

   % find faces belonging to coarseBlocks
   act1 = false([size(CG.cells.faces,1), 1]);
   act1(CS.activeCellFaces) = true;

   i       = mcolon(CG.cells.facePos(coarseBlocks    ), ...
                    CG.cells.facePos(coarseBlocks + 1) - 1) .';
   act2    = false([size(CG.cells.faces, 1), 1]);
   act2(i) = true;

   faces = unique(CG.cells.faces(act1 & act2, 1));
end
%-----------------------------------------------------------------------

function s = id(s)
s = ['updateBasisFunc:', s];
end
