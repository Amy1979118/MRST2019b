function [ fluid ] = addThermalFluidProps( fluid, varargin )
% Add thermal properties to an existing fluid structure
%
% SYNOPSIS:
%
%  fluid = addThermalFluidProps(fluid,'pn1', pv1, ...);
%
%  fluid = addThermalFluidProps(fluid,'pn1', pv1, 'useEOS', true, 'brine', true, ...);
%
% PARAMETERS:
%   fluid   - fluid structure created with initSimpleADIFluid (or other).
%
%   cp      - fluid heat capacity in J kg-1 K-1. typically 4.2e3. Used for
%             the evaluation of the internal energy and enthalpy.
%
%   lambdaF - fluid heat conductivity in W m-1 K-1. typically 0.6.
%
%   useEOS  - logical. By default the density/viscosity formulation
%              of Spivey is used.
%
%   brine   - logical. Used for the EOS and p,T,c dependency of properties
%             if salt is present or not.
%
%   dNaCl   - salt molecular diffusivity in m2s-1.
% 
% 
%   rho     - density at reservoir condition. Can be given as an handle 
%             function.
% 
%   useBFactor - logical. Used if we want bW definition factor instead 
%                of rhoW 
% 
%
% RETURNS:
%   fluid - updated fluid structure containing the following functions
%           and properties.
%             * bX(p,T,c)     - inverse formation volume factor
%             * muX(p,T,c)    - viscosity functions (constant)
%             * uW(p,T,c)     - internal energy
%             * hW(p,T,c)     - enthalpy
%             * lambdaF       - fluid conductivity
%             * dNaCl         - Salt molecular diffusivity
%
%
% SEE ALSO:
%
% 'initSimpleADIFluid', 'initSimpleThermalADIFluid'
Watt = joule/second;
opt = struct('cp'               , 4.2*joule/(Kelvin*gram), ...
             'lambdaF'          , 0.6*Watt/(meter*Kelvin), ...
             'useEOS'           , false                  , ...
             'brine'            , false                  , ...
             'dNaCl'            , 1e-6                   , ...
             'rho'              , fluid.rhoWS            , ...
             'useBFactor'       , false                  );

opt = merge_options(opt, varargin{:});

if isa(opt.rho, 'function_handle')
    asser(~opt.useBFactor, 'You cannnot use b factors and provide rho as a function handle');
end

fluid.lambdaF = opt.lambdaF;

fNames = fieldnames(fluid);
phases = '';
for fNo = 1:numel(fNames)
    fn = fNames{fNo};
    if strcmpi(fn(1:2), 'mu')
        phases = [phases, fn(3)];
    end
end
nPh   = numel(phases);

fluid.dNaCl = opt.dNaCl;

if any(strcmpi(phases, 'W'))
    if opt.useBFactor
        if ~opt.brine
            fluid.rhoW = @(p,T) fluid.bW(p).*fluid.rhoWS;
        else
            fluid.rhoW = @(p,T,c) fluid.bW(p).*fluid.rhoWS;
        end
    else
        if ~opt.brine
            fluid.rhoW = @(p,T) opt.rho;
        else 
            fluid.rhoW = @(p,T,c) opt.rho;
        end 
    end
end

if any(strcmpi(phases, 'W')) && opt.useEOS
    if ~opt.brine
        fluid.rhoW  = @(p,T) density_pure_water(p,T);
        fluid.muW   = @(p,T) viscosity_pure_water(p,T);
    else
        fluid.rhoW  = @(p,T,c) density_brine(p,T,c);
        fluid.muW   = @(p,T,c) viscosity_brine(p,T,c);
    end
end

if ~opt.useBFactor
    if ~opt.brine
        fluid.rhoWS = fluid.rhoW(1*atm, 273.15 + 20);
    else
        fluid.rhoWS = fluid.rhoW(1*atm, 273.15 + 20, 0);
    end
end

names = upper(phases);
cp    = opt.cp;
for phNo = 1:nPh
    n = names(phNo);
    if ~opt.brine
        fluid.(['h', n]) = @(p,T) cp(phNo)*T + p./feval(fluid.(['rho', n]), p,T);
        fluid.(['u', n]) = @(p,T) cp(phNo)*T;
    else
        fluid.(['h', n]) = @(p,T,c) cp(phNo)*T + p./feval(fluid.(['rho', n]), p,T,c);
        fluid.(['u', n]) = @(p,T,c) cp(phNo)*T;
    end
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
