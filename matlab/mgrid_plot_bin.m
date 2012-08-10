function [time, h_w, h_fill] = mgrid_plot_bin(n_t, n_grid, w_level, ...
    plot_body, h_fill, line_color)
% Plots straight from a binary file (i.e. ./output/ibfsxxxxxxx.var)

[x, y, u, v, w, f_b] = mgrid_read_bin(n_t);

x = x(:,:,n_grid);
y = y(:,:,n_grid);
u = u(:,:,n_grid);
v = v(:,:,n_grid);
w = w(:,:,n_grid);

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
    h_fill = true; % Fill in contours
end

if nargin < 6
    line_color = ''; % Contour line colors
end

h_fill_body = false;

%% INPUTS

poly_color = [0.5 0.5 0.5]; % Gray for the polygon

%% READ INPUT FILE AND GET PARAMETERS
s = mgrid_read_input('./input/ibfs.inp');

% Some things needed for read2d_DB1()
nx = s.M - 1;
ny = s.N - 1;
time = n_t*s.DT;
n_start = s.ISTART;

% Get the min and max and make a linear space for w
if auto == true    
    w_level = linspace(min(min(w)), max(max(w)), n_levels); %#ok<NASGU>
    w_level = level_scale_factor*w_level;
end

% Get axis limits
a_lim = [min(min(x)), max(max(x)), min(min(y)), max(max(y))];

%% PLOT

%Plot the vorticity
h_fig = gcf;
set(h_fig, 'Position', [0 1200 1200 600]);
set(h_fig, 'Color', [1 1 1]);
[cons, h_w] = contour(x, y, w, w_level, line_color);

if h_fill == true
    colorbar;
    % Here is a nice fix to fill the last contour on a filled plot
    set(gca, 'Color', [0 0 .5]);
    set(h_w, 'Fill', 'on');
end
axis equal; 
axis(a_lim);
set(gcf, 'InvertHardCopy', 'off');
% title_string = sprintf('$t = %4.4f$', time);
% title(title_string)

if plot_body == true
    % Plot the body
    [x_b, y_b, n_body, n_points] = mgrid_read_geom('./input/geom.inp');
    hold on;
    if h_fill_body == true        
        h_fill = fill(x_b, y_b, poly_color); % Fill the body!
        set(h_fill, 'EdgeColor', poly_color);
    else
        plot(x_b, y_b, 'k-o', 'LineWidth', 2);
    end
    hold off;
end

display('Done plotting...')

