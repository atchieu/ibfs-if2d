function mgrid_animate_vor_bin()
%MGRID_ANIMATE_VOR Animate MGRID simulation vorticty plot.
%   MGRID_ANIMATE_VOR() animates the vorticity that is associated with the
%   corresponding input file in the current directory. It will create a
%   output directory called ./graphics and place all the post processed
%   laprint files in that directory. Once done, one must invoke to the
%   batchtextopdf script to convert the files to PDFs. For more details see
%   MGRID_PLOT_VOR.
%
%   See also MGRID_PLOT_VOR

read_0 = true;
plot_body = true;

s = mgrid_read_input('./input/ibfs.inp');

n_read = s.ISAVE; % How many to skip before plotting another
n_total =  s.ISTOP - s.ISTART;
i_end = n_total/n_read-1;
n_grid = 2;
w_level = -1.5:0.2:1.5;

% Make a directory called output to store the figures
!mkdir ./graphics

for i = 1:i_end+1
    n_t = i*n_read + s.ISTART;
    n_t_pad = pad(n_t, 7);
    fprintf('\nPlotting frame_t%s!!\n', n_t_pad);
    time = mgrid_plot_bin(n_t, n_grid, w_level);
    if plot_body == true
        hold on;
%         mgrid_plot_body(time)
%         mgrid_plot_arrow(time)
        hold off;
    end
    title_string = sprintf('$t = %4.2f$', time);
    title(title_string);
    output_name = sprintf('./graphics/frame_t%s', n_t_pad);
    laprint(1, output_name, 'width', 10.0);
end


function mgrid_plot_body(time)
% Just an auxilary function to plot the body, you make this up as you go
x = [-0.5 0.5]';
y = [ 0.0 0.0]';

poly_color = [0.5 0.5 0.5]; % Gray for the polygon

theta = 0:0.01:2*pi;
x = real(0.5*exp(i*theta));
y = imag(0.5*exp(i*theta));

theta = time;
x_p =  x*cos(theta);
y_p =  x*sin(theta);

% plot(x, y, 'k', 'LineWidth', 1);
h_fill = fill(x, y, poly_color); % Fill the body!
set(h_fill, 'EdgeColor', poly_color);


function mgrid_plot_arrow(time)

theta = 2./pi * sin(pi*time);

x = [-1.5, -1.0, -1.15, -1.15, -1.0];
y = [ 0.0,  0.0, 0.03, -0.03, 0.0];

x_p =  x*cos(theta) + y*sin(theta);
y_p = -x*sin(theta) + y*cos(theta);

% plot(x_p, y_p, k);
h_fill = fill(x_p, y_p, [0 0 0]); % Fill the body!
set(h_fill, 'EdgeColor', [0 0 0]);