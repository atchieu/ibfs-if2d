%% ------------------------------------------------------------------------
% function read_geom(file)
%
% Reads the Lagrangian points provided by the appropriate file.
%
%   file        Input file name.
%--------------------------------------------------------------------------
function [x, y, n_body, n_points] = mgrid_read_geom(file)

% Check the correct input arguments
if nargin < 1
  error('Function requires one input argument');
elseif ~ischar(file)
  error('Input must be a string representing a filename');
end

% Open file
fid = fopen(file);
if fid == -1
    error('Could not find the file %s', file)
end

n_points = fscanf(fid, '%d', 1);
data = fscanf(fid, '%f %f %d', [3, inf]);
data = data';
x = data(:, 1);
y = data(:, 2);
n_body = data(:, 3);

fclose(fid);

display('Done reading geometry file...')