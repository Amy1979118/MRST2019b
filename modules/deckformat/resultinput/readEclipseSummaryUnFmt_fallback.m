function [smry, smspec] = readEclipseSummaryUnFmt(prefix)
%Read unformatted (binary) ECLIPSE summary data
%
% SYNOPSIS:
%   [summary, smspec] = readEclipseSummaryUnFmt(prefix)
%
% PARAMETERS:
%   prefix - Path-name prefix from which to construct list of summary file
%            names.  Specifically, this function reads files which match
%            the regular expressions
%
%                [prefix, '\.SMSPEC']  (unformatted summary specification)
%                [prefix, '\.S\d{4}']  (unformatted summary files)
%
%            Use function 'readEclipseSummaryFmt' to read formatted
%            (text/ASCII) summary data.
%
% RETURNS:
%   summary - Summary data structure.  All MINISTEP results of an ECLIPSE
%             run concatenated together.  Also includes an additional
%             sub-structure, UNITS, whose fields contain a textual
%             representation of unit of measurement of the corresponding
%             summary field such as
%
%                 summary.UNITS.TIME  = 'DAYS'
%                 summary.UNITS.WBHP  = 'BARSA'
%                 summary.UNITS.WOPR  = 'SM3/DAY'
%                 summary.UNITS.YEARS = 'YEARS'
%                 summary.UNITS.TCPU  = 'SECONDS'
%
%   smspec  - Summary specifiction obtained from the '.SMSPEC' file.
%             Contains MINISTEP times and all summary vectors declared in
%             the run deck.  Additionally contains a field, '.RptTime' that
%             specifies the times at which restart files (3D data) is
%             reported.
%
% SEE ALSO:
%   `readEclipseSummaryFmt`.

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


   is_open_pre = fopen('all');

   [dname, fp] = fileparts(prefix);
   if isempty(dname),
      dname = '.';
   end

   smspec = readEclipseOutputFileUnFmt([prefix, '.SMSPEC']);
   summaryfiles = matchResultFiles(dname, [fp, '\.S\d{4}']);

   smry = readEclipseSummary(summaryfiles, smspec,            ...
                             @(fn) fopen(fn, 'r', 'ieee-be'), ...
                             @readFieldUnFmt);

   is_open_post = fopen('all');

   assert (all(size(is_open_pre) == size(is_open_post)) && ...
           all(is_open_pre == is_open_post),               ...
           'Summary reader leaks file identifiers.');
end
