function [ eqs ] = addEnthalpyContributionGeothermal(model, state, eqs, names, src,p, T, c, drivingForces)
% Add the enthalpy contribution for boundaries and source terms to equations
%
% SYNOPSIS:
%   [eqs] = ...
%    addEnthalpyContributionGeothermal(model, ....
%                                         eqs, names, src, p, T, c,
%                                                drivingForces )
%
% PARAMETERS:
%   model  - Class instance.
%
%   eqs    - Cell array of equations that are to be updated.
%
%   names  - The names of the equations to be updated. If
%            phase-pseudocomponents are to be used, the names must
%            correspond to some combination of "water", "NaCl", "temperature"
%            if no special component treatment is to be introduced.
%
%   types  - Cell array with the types of "eqs". Note that these
%            types must be 'cell' where source terms is to be added.
%
%   src    - Struct containing all the different source terms that
%            were computed and added to the equations.
%
%   p      - Cell array of phase pressures.
%
%   T      - Cell array of phase temperature.
%
%   c      - Cell array of phase NaCl mass fraction.
%
%   forces - DrivingForces struct (see 'getValidDrivingForces')
%            containing (possibily empty) 'src' and 'bc' fields.
%
% RETURNS:
%   eqs   - Equations with corresponding source terms and BCs added.


% check that the mass flux is computed
assert(~isempty(src.bc.phaseMass{1}),'the mass flux is not computed. Pressure boundary are not imposed')

% get all needed variables
G       =  model.G;
op      =  model.operators;
f       =  model.fluid;
bc      =  drivingForces.bc;
bcT     =  drivingForces.bcT;
Tr_bcd  =  op.T_all(bc.face);

T_heat_R_all = op.T_heat_R_all;
T_heat_F_all = op.T_heat_F_all;
T_heat_all   = T_heat_F_all + T_heat_R_all;


% Get water density and mobility
[rho, mob]   = model.getProps(state, 'Density', 'Mobility');
rhoW  = rho{1}; 
mobW  = mob{1};        

id_eq   =  strcmpi(names,'temperature');
bcCells =  src.bc.sourceCells;

if(~isempty(bcCells))
    %  get the value of P at the faces of the pressure boundary (will be needed to compute the enthalpy)
    %  the flux for pressure is defined in the function
    %  "getBoundaryConditionFluxesAD" as : q = T.*mob.*dP;
    %  with dP = Pbcf - Pbcc + rhoBC .*dzbc; T: trans. and mob = 1/muW.
    %
    
    PBCFaces = bc.face;
    % gravity gradient per bc face
    if any(strcmpi(G.type,'topSurfaceGrid'))
        dzbc = model.gravity(3) * (G.cells.z(bcCells) - G.faces.z(bc.face));
    else
        g    =  model.getGravityVector();
        dz   =  G.cells.centroids(bcCells,:) - G.faces.centroids(PBCFaces,:);
        dzbc =  dz*g';
    end
    
    isP = reshape(strcmpi(bc.type, 'pressure'),[],1);
    isF = reshape(strcmpi(bc.type, 'flux'),[],1);
    
    pbcf = p(bcCells);
    if isprop(model, 'NaCl') 
        cbcf = c(bcCells);
    end 
    
    if any(isP)
        pbcf(isP) = bc.value(isP);
    end
    
    if any(isF)
        mobWbc    =  mobW(bcCells);
        rhoWbc    =  rhoW(bcCells);
        dP        =  bc.value(isF)./mobWbc(isF)./Tr_bcd(isF);
        pbcf(isF) =  dP + pbcf(isF) - rhoWbc(isF).*dzbc(isF);
    end
    
    % Get the value for temperature at BCs (initialisation)
    % by default take the value at the cell centroid
    tbcf   =  T(bcCells);
    
    % Sanity check : if any bcT are prescribed
    if(~isempty(bcT))
        Tr_bch  =  T_heat_all(bcT.face);
        % check if we have info on temp at BCs and update tbcf
        TBCFaces    =  bcT.face;
        ind         =  find(ismember(PBCFaces, TBCFaces));
        
        % update enthalpy structure based on the temperature condition (fixed or flux)
        if any(ind)
            faces =  PBCFaces(ind);
            idft  =  find(ismember(TBCFaces, faces));
            isT   =  reshape(strcmpi(bcT.type(idft),'temperature'), [],1);
            isTF  =  reshape(strcmpi(bcT.type(idft),'Tflux'), [],1);
            
            inflow = value(src.bc.phaseMass{1}) > 0;
            
            if any(isT) % there is a 'temperature' condition
                idxt       =  idft(isT);                                     % mapping in Tbc vector, indexes that have an info on P.
                tfaces     =  bcT.face(idxt);                                % which faces have an info on P and a fixed 'temperature'
                idxp       =  ismember(PBCFaces, tfaces);                    % get the indexes of these faces in the vector PBCFaces
                tbcf(idxp) =  inflow(isT).*bcT.value(idxt) + ~inflow(isT).*tbcf(isT);                               % update the temperature face at BC
            end
            
            if any(isTF) % there is a 'heat flux' condition
                idxt       =  idft(isTF);                                    % mapping in Tbc vector, indexes that have an info on P
                tfaces     =  bcT.face(idxt);                                % which faces have an info on P and a fixed 'temperature'
                idxp       =  ismember(PBCFaces, tfaces);                    % get the indexes of these faces in the vector PBCFaces
                dT         =  bcT.value(idxt)./Tr_bch(idxt);                 % update the temperature face at BC
                tbcf(idxp) =  dT + tbcf(idxp);
            end
            
        end
    end
    
    if isprop(model, 'NaCl') 
        hWbcf  =  f.hW(pbcf, tbcf, cbcf);
    else 
        hWbcf  =  f.hW(pbcf, tbcf);
    end 
     
    
    qBC         =   src.bc.phaseMass{1};
    qtbc        =   src.bc.mapping*(hWbcf.*qBC);
    
    eqs{id_eq}(bcCells) = eqs{id_eq}(bcCells) - qtbc;
end
end

%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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
