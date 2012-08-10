function mgrid_animate_vor()
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
plot_body = false;

s = mgrid_read_input('./input/ibfs.inp');

n_read = s.ISAVE; % How many to skip before plotting another
n_total =  s.ISTOP - s.ISTART;
i_end = n_total/n_read-1;
n_grid = 4;
w_level = 0.05:.05:1.4;

% Make a directory called output to store the figures
!mkdir ./graphics

for i = 1:i_end+1
    n_t = i*n_read;
    n_t_pad = pad(n_t, 7);
    fprintf('\nPlotting frame_t%s!!\n', n_t_pad);
    time = mgrid_plot_vor(n_t, n_grid, w_level);
    if plot_body == true
        hold on;
        mgrid_plot_body(time)
        hold off;
    end
    output_name = sprintf('./graphics/frame_t%s', n_t_pad);
    laprint(1, output_name, 'width', 10.0);
end


function mgrid_plot_body(time)
% Just an auxilary function to plot the body, you make this up as you go
x = [-0.5 0.5]';
y = [ 0.0 0.0]';

time = 0;

theta = pi/5 * time;
x_p =  x*cos(theta);
y_p =  x*sin(theta);

plot(x_p, y_p, 'k', 'LineWidth', 2);
