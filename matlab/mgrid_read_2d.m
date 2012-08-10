function [x, y, u, v, w, s] = mgrid_read_2d(n_t, n_grid)

%% READ INPUT FILE AND GET PARAMETERS
s = mgrid_read_input('./input/ibfs.inp');

% Some things needed for read2d_DB1()
nx = s.M - 1;
ny = s.N - 1;

% Create a readable filename
n_t_pad = pad(n_t, 7); % create a padded version of n_t
filename = sprintf('./post/ibfs_g%dt%s_vel.dat', n_grid, n_t_pad);

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

display('Done reading .dat file...')