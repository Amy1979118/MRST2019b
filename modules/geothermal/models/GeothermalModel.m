classdef GeothermalModel < WaterModel 
    %
    % DESCRIPTION:
    %   General model for geothermal problems. Contains BCs and caping for
    %   the heat equation. GeothermalModel is a subclass of the WaterModel.
    %
    % REQUIRED PARAMETERS:
    %   G     - Simulation grid.
    %
    %   rock  - Valid rock used for the model.
    %
    %   fluid - Fluid model used for the model.
    %
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to the specified value.
    %
    % RETURNS:
    %   Class instance.
    %
    % SEE ALSO:
    %   `GeothermalWaterModel`, `GeothermalWaterNaClModel`
    
    properties
        thermal = true;
        maximumTemperature = Inf;
        minimumTemperature = -Inf;
    end
    
    methods
        function model = GeothermalModel(G, rock, fluid, varargin)
            model       =   model@WaterModel(G, rock, fluid);
            model       =   merge_options(model, varargin{:});
            lambdaR     =   rock.lambdaR; 
            lambdaF     =   fluid.lambdaF; 
            poro        =   rock.poro;

            % transmissibilities for fluid and rock            
            LR = (1-poro).*lambdaR; 
            LF = poro.*lambdaF;              
            rock_heat  = struct('perm',LR); 
            fluid_heat = struct('perm',LF);            
            T_heat_R   = computeTrans(G,rock_heat);
            T_heat_F   = computeTrans(G,fluid_heat);
    
            % harmonic average of transmissibilities
            cf          =   G.cells.faces(:,1);
            nf          =   G.faces.num;
            T_heat_R      =   1 ./ accumarray(cf, 1./T_heat_R, [nf, 1]);
            T_heat_F      =   1 ./ accumarray(cf, 1./T_heat_F, [nf, 1]);
            model.operators.T_heat_R_all      =   T_heat_R;
            model.operators.T_heat_F_all      =   T_heat_F;
            
            intInx                            =   all(G.faces.neighbors ~= 0, 2);
            model.operators.T_heat_R          =   T_heat_R(intInx);
            model.operators.T_heat_F          =   T_heat_F(intInx);
            
            % add thermal state functions for property evaluation
            model.FlowPropertyFunctions = ThermalFlowPropertyFunctions(model);
            model.FluxDiscretization    = ThermalFluxDiscretization(model);
            
            
        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@ReservoirModel(model);
            forces.bcT = [];
            forces.srcT = [];
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state = model.updateStateFromIncrement(state, dx, problem, 'temperature', 0.2, inf);
            state = model.capProperty(state, 'temperature', model.minimumTemperature, model.maximumTemperature);
            [problem.primaryVariables, removed] = model.stripVars(problem.primaryVariables, 'temperature');
            dx = dx(~removed);

            [state, report] = updateState@WaterModel(model, state, problem, dx, drivingForces);
        end
        
        function rhoS = getSurfaceDensities(model)
            active = model.getActivePhases();
            props = {'rhoWS', 'rhoOS', 'rhoGS'};
            rhoS = cellfun(@(x) model.fluid.(x), props(active));
        end
        
         function [ eqs ] = addThermalBoundaryConditionsGeothermal( model, eqs, names, p, T, drivingForces )
           
        % Add the thermal boundary conditions to equations
        %
        % SYNOPSIS:
        %   [eqs] = .....
        %    addThermalBoundaryConditionsGeothermal(model, ....
        %                                           eqs, names,  p, T,
        %                                           drivingForces )
        %
        % PARAMETERS:
        %   model  - Class instance.
        %   eqs    - Cell array of equations that are to be updated.
        %
        %   names  - The names of the equations to be updated. If
        %            phase-pseudocomponents are to be used, the names must
        %            correspond to some combination of "water", 
        %            "temperature", "NaCl" if no special component 
        %            treatment is to be introduced.
        %
        %   p      - Cell array of phase pressures.
        %
        %   T      - Cell array of phase temperature.
        %
        %
        %   forces - DrivingForces struct (see `getValidDrivingForces`)
        %            containing (possibily empty) `src` and `bc` fields.
        %
        % RETURNS:
        %   eqs   - Equations with corresponding source terms and BCs added.
            
            
            
            hasBC   = isfield(drivingForces, 'bcT')  && ~isempty(drivingForces.bcT);
            id_eq = strcmpi(names,'temperature');
            
            if hasBC
                 
                 bcT        = drivingForces.bcT;
                 G          = model.G;
                 N          = G.faces.neighbors(bcT.face,:);
                 BCcells    = sum(N,2);
                 op         = model.operators; 

                 T_heat_R_all  = op.T_heat_R_all;
                 T_heat_F_all  = op.T_heat_F_all;
                 T_heat_all    = T_heat_F_all + T_heat_R_all;
                 T_h           = T_heat_all(bcT.face);
              
                 % read the type of BCs: temperature, flux and enthalpy
                 % temperature and flux BCs are conditions for the rock
                 % enthalpy BC is a condition for the fluids (if we want to set the fluid
                 % with a different temperature
                 isT     =   reshape(strcmpi(bcT.type,'temperature'),[],1);
                 isTF    =   reshape(strcmpi(bcT.type,'Tflux'),[],1);
                 
                 q_s = zeros(numel(BCcells), 1);
                 q_s = model.AutoDiffBackend.convertToAD(q_s, T);
                
                % Treat temperature BCs
                % Get the contribution from the rock
                % heat flux qbT for the different type of BCs:
                % Dirichlet BC
                if any(isT)
                    subs        =    isT;
                    bcc_id      =    BCcells(subs);
                    Tbc         =    T(bcc_id);
                    Tbf         =    bcT.value(subs);
                    T_h_id      =    T_h(subs);                     
                    dT          =    Tbf - Tbc;
                    qbT         =    T_h_id.*(dT);
                    
                    q_s(subs) =    qbT;
                end
                
                % Neumann BC
                if any(isTF)
                    subs      =    isTF;
                    qbT       =    bcT.value(subs);
                    q_s(subs) =    qbT;
                end
                S = sparse(BCcells, (1:numel(bcT.face))', 1, G.cells.num, numel(bcT.face));
                
                % Update equation:
                eqs{id_eq} = eqs{id_eq} - S*q_s;
             
            end
         end
     
         function[ eqs, srcT ] = addThermalSourcesGeothermal(model, eqs, names, T, drivingForces ) %#ok
        % Add in the thermal source terms to equations 
        % 
        % SYNOPSIS:
        %   [eqs] = .....
        %    addThermalSourcesGeothermal(model, ....
        %                                  eqs, names,  p, T,
        %                                  drivingForces )
        %
        % PARAMETERS:
        %   model  - Class instance.
        %   eqs    - Cell array of equations that are to be updated.
        %
        %   names  - The names of the equations to be updated. If
        %            phase-pseudocomponents are to be used, the names must
        %            correspond to some combination of "water", 
        %            "temperature", "NaCl" if no special component 
        %            treatment is to be introduced.
        %
        %   T      - Cell array of phase temperature.
        % 
        %   forces - DrivingForces struct (see `getValidDrivingForces`)
        %            containing (possibily empty) `src` and `bc` fields.
        %
        % RETURNS:
        %   eqs   - Equations with corresponding source terms added.
       
            hasSRC   = isfield(drivingForces, 'srcT')  && ~isempty(drivingForces.srcT);
            id_eq = strcmpi(names,'temperature');
            
            srcT = []; 
            if hasSRC
                   srcT         = drivingForces.srcT; 
                   SRCcells     = srcT.cell;
                   Tbcell       = T(SRCcells); 
                  
                   eqs{id_eq}(SRCcells) = Tbcell - srcT.value;                    
            end
         end 
         
              
       function [eqs, names, types, wellSol] = insertWellEquations(model, eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, T,components, dt, opt)
        % Utility function for setting up the well equations and adding
        % source terms for black-oil like models. Note that this currently
        % assumes that the first nPh equations are the conservation
        % equations, according to the canonical MRST W-O-G ordering,
        fm = model.FacilityModel;
        active=model.getActivePhases;
        nPh = nnz(active);
        %assert(numel(eqs) == nPh+2);
        dissolved =[];
        [eqs, names, types, wellSol, src] = insertWellEquations@ReservoirModel(model, eqs, names, ...
            types, wellSol0, wellSol, ...
            wellVars, ...
            wellMap, p, mob, rho, ...
            dissolved, components, ...
            dt, opt);


        actWellIx = fm.getIndicesOfActiveWells(wellSol); 
        nw = numel(actWellIx);

        % sanity check (if no well structure)
            if (~isempty(src))

                % update water enthalpy (hW) for the well equation
                % compute hW with the pressure and temperature info
                % if a Tinj is prescribed use it to compute hW

                if (~isempty(src.phaseMass))
                    f   = model.fluid;
                    idt = strcmpi(names,'temperature');
                    src2well = nan(model.G.cells.num,1);
                    src2well(src.sourceCells) = 1:numel(src.sourceCells);

                    for i = 1:nw
                        wellNo = actWellIx(i);
                        wm     = fm.WellModels{wellNo};
                        ws     = wellSol(wellNo);
                        wc     = wm.W.cells;
                        tw     = T(wc);
                        cqsw   = src.phaseMass{1}(src2well(wc));

                        tw(cqsw>0) = wm.W.T;
                        pw     = p(wc);
                        hw     = f.hW(pw,tw);
                        eqs{idt}(wc) = eqs{idt}(wc) - hw.*cqsw;

                    end
                end
            end
       end
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