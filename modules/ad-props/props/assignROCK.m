function f = assignROCK(f, rock, reg)
    ntpvt = reg.pvt;
    cR   = rock(:, 2)';
    pRef = rock(:, 1)';
    if ntpvt == 1
        f.pvMultR = @(p) pvMult(p, cR, pRef);
    else
        f.pvMultR = cell(1, ntpvt);
        for i = 1:ntpvt
            f.pvMultR{i} = @(p) pvMult(p, cR(i), pRef(i));
        end
    end
end

function v = pvMult(p, cR, pRef)
    x = cR.*(p-pRef);
    v = 1 + x + 0.5*(x.^2);
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
