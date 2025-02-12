function [info, present] = dataset_npd_bjarmeland()
% Info function for Frigg dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

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
    [info, present] = datasetInfoStruct(...
        'name', 'NPD_Bjarmeland', ...
        'website', 'https://www.npd.no/fakta/co2-handtering-ccs/lagringsmodeller/', ...
        'fileurl', 'https://www.npd.no/globalassets/1-npd/fakta/co-to/lagringsmodeller/barentshavet/bjarmeland.zip', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', true, ...
        'cells',   nan, ...
        'examples', {}, ...
        'description', 'Barents sea CO2 storage model from the Norwegian Petroleum Directorate (NPD)',...
        'filesize',    370, ...
        'modelType', 'Eclipse deck (CO2STORE).' ...
         );
end
