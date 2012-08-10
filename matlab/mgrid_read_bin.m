function [x, y, u, v, w, f_b] = mgrid_read_bin(n_t)
% Read straight from a binary file (i.e. ./output/ibfsxxxxxxx.var)

%% READ INPUT FILE, GEOMETRY FILE, AND GET PARAMETERS
s = mgrid_read_input('./input/ibfs.inp');
[x_b, y_b, n_body, n_f] = mgrid_read_geom('./input/geom.inp'); 

m = s.M;
n = s.N;
n_k = s.MGRIDLEV;
delta = s.LEN/s.M;

% Check the maximum MGRIDLEV
n_levmax = floor( log(min(m,n)/2)/log(2.) ) - 1;
n_k = min(s.MGRIDLEV, n_levmax);

% We need to read omega, fb, rhs_old, q, and q0
%
%   omega and rhs_old is size (2:m, 2:n, n_k)
%   fb is size (n_f)
%   q and q0 is size (n_q, n_k)
n_q = (m+1)*n + (n+1)*m;
n_f = 2*n_f; % For both x and y
n_w = (m-1)*(n-1);

%% OPEN FILE AND READ A RIDICULOUS AMOUNT OF INFORMATION

% Create a readable filename
n_t_pad = pad(n_t, 7); % create a padded version of n_t
filename = sprintf('./output/ibfs%s.var', n_t_pad);

fid = fopen(filename, 'r');
discard = fread(fid, 1, 'int32'); % Must discard (?), FORTRAN problem
w_col = fread(fid, n_w*n_k, 'float64');
f_b = fread(fid, n_f, 'float64');
discard = fread(fid, n_w*n_k, 'float64');
q_col = fread(fid, n_q*n_k, 'float64');
q0_col = fread(fid, n_q*n_k, 'float64');

w = reshape(w_col, m-1, n-1, n_k);
q = reshape(q_col, n_q, n_k);
q0 = reshape(q0_col, n_q, n_k);

fclose(fid);

display('Done reading binary file...')

%% CREATE GRID

% Poor man's way of keeping indices (otherwise rewrite all these column
% major routines in loops)
u_i = zeros(m+1,n);
v_i = zeros(m,n+1);

next = 0;
for j = 1:n
    for i = 1:m+1
        next = next + 1;
        u_i(i,j) = next;
    end
end

for j=1:n+1
    for i = 1:m
        next = next + 1;
        v_i(i,j) = next;
    end
end

% Allocate space
u = zeros(m,n,n_k);
v = zeros(m,n,n_k);
x = zeros(m,n,n_k);
y = zeros(m,n,n_k);

for k = 1:n_k
    % Multi-domain factors and offsets
    fac = 2^(k-1);
    off_x = (fac-1.)/2. * s.LEN;
    off_y = (fac-1.)/2. * delta*n;
    del = fac*delta;

    % Scale vorticity
    w(:,:,k) = w(:,:,k)./del^2;

    for j = 2:n
        for i = 2:m
            u(i,j,k) = 0.5*( q(u_i(i,j),k)  + q(u_i(i,j-1),k) + ...
                             q0(u_i(i,j),k) + q0(u_i(i,j-1),k) ) / del;
            v(i,j,k) = 0.5*( q(v_i(i,j),k)  + q(v_i(i-1,j),k) + ...
                             q0(v_i(i,j),k) + q0(v_i(i-1,j),k) ) / del;
            x(i,j,k) = del*(i-1) - off_x - s.OFFSETX;
            y(i,j,k) = del*(j-1) - off_y - s.OFFSETY;
        end
    end
end

u = u(2:m,2:n,:);
v = v(2:m,2:n,:);
x = x(2:m,2:n,:);
y = y(2:m,2:n,:);

n_total = 2*n_w*n_k + n_f + 2*n_q*n_k;
