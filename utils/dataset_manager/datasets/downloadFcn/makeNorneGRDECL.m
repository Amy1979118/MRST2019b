function ok = makeNorneGRDECL
%Create containing datafile for subset of Norne simulation model
%
% SYNOPSIS:
%   ok = makeNorneGRDECL
%
% DESCRIPTION:
%   This function ensures existence of the file
%
%       fullfile(getDatasetPath('norne'), 'NORNE.GRDECL')
%
%   that contains INCLUDE statements for an appropriate subset of the Norne
%   simulation model.
%
% PARAMETERS:
%   None.  This function operates on constant data.
%
% RETURNS:
%   ok - Whether or not we were able to create the datafile.
%
% SEE ALSO:
%   `getDatasetPath`.

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

   ok = grdecl_available();

   if (~ ok) && makeNorneSubsetAvailable(),
      ok = create_grdecl();
   end
end

%--------------------------------------------------------------------------

function ok = grdecl_available()
   ok = exist(grdecl_filename(), 'file') == 2;
end

%--------------------------------------------------------------------------

function ok = create_grdecl()
   assert (~ grdecl_available(), 'Internal Error');

   [fid, msg] = fopen(grdecl_filename(), 'wt');

   if fid < 0,
      warning('Open:Fail', 'Failed to open GRDECL output file: %s\n', msg);

      ok = false;

      return;
   end

   write_grdecl(fid);

   fclose(fid);

   ok = true;
end

%--------------------------------------------------------------------------

function file = grdecl_filename()
   ndir = getDatasetPath('norne', 'skipAvailableCheck', true);
   file = fullfile(ndir, 'NORNE.GRDECL');
end

%--------------------------------------------------------------------------

function write_grdecl(fid)
   fprintf(fid, [ ...
'INCLUDE\n', ...
'  ''INCLUDE/IRAP_1005.GRDECL'' /\n\n', ...
'INCLUDE\n', ...
'  ''INCLUDE/ACTNUM_0704.prop'' /\n\n', ...
'INCLUDE\n', ...
'  ''INCLUDE/PERM_0704.prop'' /\n\n', ...
'INCLUDE\n', ...
'  ''INCLUDE/PORO_0704.prop'' /\n\n', ...
'INCLUDE\n', ...
'  ''INCLUDE/NTG_0704.prop'' /\n\n', ...
'COPY\n', ...
'   ''PERMX'' ''PERMY'' /\n', ...
'   ''PERMX'' ''PERMZ'' /\n', ...
'/\n\n', ...
'MULTIPLY\n', ...
'   ''PERMZ'' 0.2    1 46 1 112  1  1 /\n', ...
'   ''PERMZ'' 0.04   1 46 1 112  2  2 /\n', ...
'   ''PERMZ'' 0.25   1 46 1 112  3  3 /\n', ...
'   ''PERMZ'' 0.0    1 46 1 112  4  4 /\n', ...
'   ''PERMZ'' 0.13   1 46 1 112  5  5 /\n', ...
'   ''PERMZ'' 0.13   1 46 1 112  6  6 /\n', ...
'   ''PERMZ'' 0.13   1 46 1 112  7  7 /\n', ...
'   ''PERMZ'' 0.13   1 46 1 112  8  8 /\n', ...
'   ''PERMZ'' 0.09   1 46 1 112  9  9 /\n', ...
'   ''PERMZ'' 0.07   1 46 1 112 10 10 /\n', ...
'   ''PERMZ'' 0.19   1 46 1 112 11 11 /\n', ...
'   ''PERMZ'' 0.13   1 46 1 112 12 12 /\n', ...
'   ''PERMZ'' 0.64   1 46 1 112 13 13 /\n', ...
'   ''PERMZ'' 0.64   1 46 1 112 14 14 /\n', ...
'   ''PERMZ'' 0.64   1 46 1 112 15 15 /\n', ...
'   ''PERMZ'' 0.64   1 46 1 112 16 16 /\n', ...
'   ''PERMZ'' 0.64   1 46 1 112 17 17 /\n', ...
'   ''PERMZ'' 0.016  1 46 1 112 18 18 /\n', ...
'   ''PERMZ'' 0.004  1 46 1 112 19 19 /\n', ...
'   ''PERMZ'' 0.004  1 46 1 112 20 20 /\n', ...
'   ''PERMZ'' 1.0    1 46 1 112 21 21 /\n', ...
'   ''PERMZ'' 1.0    1 46 1 112 22 22 /\n', ...
'/\n']);

   idir = fullfile(fileparts(fopen(fid)), 'INCLUDE');

   isfile = @(f) exist(fullfile(idir, f), 'file') == 2;

   if isfile('MULTZ_HM_1.INC'),
      fprintf(fid, ...
['\nINCLUDE\n  ''INCLUDE/MULTZ_HM_1.INC'' /\n\n', ...
'EQUALS\n', ...
'  ''MULTZ''    1.0      1  46  1 112   1   1  /\n', ...
'  ''MULTZ''    0.05     1  46  1 112  15  15  /\n', ...
'  ''MULTZ''    0.001    1  46  1 112  18  18  /\n', ...
'  ''MULTZ''    0.00001  1  46  1 112  20  20  /\n', ...
'/\n']);

      if isfile('MULTZ_JUN_05_MOD.INC'),
         fprintf(fid, ...
'\nINCLUDE\n  ''INCLUDE/MULTZ_JUN_05_MOD.INC'' /\n');
      end
   end
end
