function [uv, p, omega, grid] = mgrid_read_2d_mod(n_t, n_grid)

% Like mgrid_read_2d but keeps u and v together.

%% READ INPUT FILE AND GET PARAMETERS
s = mgrid_read_input('./input/ibfs.inp');

% Some things needed for read2d_DB1()
nx = s.M - 1;
ny = s.N - 1;
time = n_t*s.DT;

% Create a readable filename
n_t_pad = pad(n_t, 7); % create a padded version of n_t
filename = sprintf('./post/ibfs_g%dt%s_vel.dat', n_grid, n_t_pad);

%% READ DATA FILE
[uv, p, omega, grid] = read2d_DB1(filename, nx, ny, 0);

