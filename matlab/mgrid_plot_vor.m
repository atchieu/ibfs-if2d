function [time, h_w, h_fill] = mgrid_plot_vor(n_t, n_grid, w_level, ...
    plot_body, fillflag, line_color)
%MGRID_PLOT_VOR plots the vorticity contours given tecplot style MGRID data
%   [h_w, h_fill] = mgrid_plot_vor(n_t, n_grid, w_level) plots the
%   vorticity contours of the time step n_t on grid level n_grid in the
%   directory ./post/ (i.e. you must be on the default directory level of
%   the program. You will also need to have the input file at the location
%   ./input/ibfs.inp. If w_level is not specified then it will plot
%   automatic contours based on the maximum and minimum levels of
%   vorticity in the field. Options include
%       
%       plot_body   if true, will plot the body given in ./input/geom.inp
%                      and ./output/com.dat
%       fill        if true, will fill the contours
%       line_color  specifies a line color for contour levels
%
%   The function returns the handles of both the contours and the filled
%   object.
%   
%   See also MGRID_ANIMATE_VOR


if nargin < 3
    auto = true; % Enable automatic contour levels
    n_levels = 50; % Number of levels in the automatic contour plot
    level_scale_factor = 1.0; % Scaling the contour levels
else
    auto = false;
end

if nargin < 4
    plot_body = false;
end

if nargin < 5
    fillflag = true; % Fill in contours
end

if nargin < 6
    line_color = ''; % Contour line colors
end


%% INPUTS

poly_color = [0.5 0.5 0.5]; % Gray for the polygon

%% READ INPUT FILE AND GET PARAMETERS
s = mgrid_read_input('./input/ibfs.inp');

% Some things needed for read2d_DB1()
nx = s.M - 1;
ny = s.N - 1;
time = n_t*s.DT;
n_start = s.ISTART;

% Create a readable filename
n_t_pad = pad(n_t, 7); % create a padded version of n_t
filename = sprintf('./post/ibfs_g%dt%s_vel.dat', n_grid, n_t_pad);

%% READ THE GEOMETRY AND PLOT BODY
[x_b, y_b, n_body, n_points] = read_geom('./input/geom.inp'); %#ok<NASGU>
display('Done reading Lagrangian body points...')

%% READ DATA FILE
[uv, p, omega, grid] = read2d_DB1(filename, nx, ny, 0);
u_col = uv(:,1); v_col = uv(:,2);
x_col = grid(:,1); y_col = grid(:,2);

% Reshape into matrices
x = reshape(x_col, nx, ny); clear x_col;
y = reshape(y_col, nx, ny); clear y_col;
u = reshape(u_col, nx, ny); clear u_col; %#ok<NASGU>
v = reshape(v_col, nx, ny); clear v_col; %#ok<NASGU>
w = reshape(omega, nx, ny); clear omega;

% Get the min and max and make a linear space for w
if auto == true    
    w_level = linspace(min(min(w)), max(max(w)), n_levels); %#ok<NASGU>
    w_level = level_scale_factor*w_level;
end

% Get axis limits
a_lim = [min(min(x)), max(max(x)), min(min(y)), max(max(y))];
% a_lim = [-2, 8, -2, 2];

display('Done reading .dat file...')

%% PLOT THE VORTICITY AND BODY

% Plot the vorticity
h_fig = figure(1);
set(h_fig, 'Position', [0 200 1200 600]);
set(h_fig, 'Color', [1 1 1]);
[cons, h_w] = contour(x, y, w, w_level, line_color);

if fillflag == true
    colorbar;
    % Here is a nice fix to fill the last contour on a filled plot
    set(gca, 'Color', [0 0 .5]);
    set(h_w, 'Fill', 'on');
end
axis equal; 
axis(a_lim);
set(gcf, 'InvertHardCopy', 'off');
title_string = sprintf('$t = %4.4f$', time);
title(title_string)

if plot_body == true
    % Plot the body
%     com = load('./output/com.dat'); % Read COM
%     x_com = com(n_t+1-n_start, 2);
%     y_com = com(n_t+1-n_start, 3);

    x_com = 0;
    y_com = 0;

    x_b = x_b + x_com;
    y_b = y_b + y_com;

    figure(1); hold on;
    h_fill = fill(x_b, y_b, poly_color); % Fill the body!
    set(h_fill, 'EdgeColor', poly_color);
    figure(1); hold off;
end

display('Done plotting...')


%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
% function read_geom(file)
%
% Reads the Lagrangian points provided by the appropriate file.
%
%   file        Input file name.
%--------------------------------------------------------------------------
function [x, y, n_body, n_points] = read_geom(file)

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
