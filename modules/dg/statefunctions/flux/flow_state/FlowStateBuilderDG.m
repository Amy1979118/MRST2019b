classdef FlowStateBuilderDG < FlowStateBuilder
    
    methods
        function flowState = build(builder, fd, model, state, state0, dt)
            
            flowState = state.faceStateDG;
            propfn = model.getStateFunctionGroupings();
            for i = 1:numel(propfn)
                p = propfn{i};
                struct_name = p.getStateFunctionContainerName();
                if isfield(state, struct_name)
                    flowState = rmfield(flowState, struct_name);
                end
            end
            flowState = model.initStateFunctionContainers(flowState);
%             state = model.validateState(state);
            flowState.type = 'face';
            if ~isfield(flowState, 'cells')
                flowState.cells = (1:model.G.cells.num)';
                flowState.faces = (1:model.G.faces.num)';
            end
        end
    end
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
