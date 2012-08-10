function mgrid_plot_force()

plot_options = 'b';

s = mgrid_read_input('./input/ibfs.inp');

fid = fopen('./output/force.dat');
forces = fscanf(fid, '%d %f %f', [3 inf])';
display('Done reading force.dat file...');

dt = s.DT;

figure(1);
subplot(2,1,1);
plot(dt*forces(:,1), forces(:,2), plot_options)
ylabel('$C_X$'); 
% set(gca, 'YLim', [1, 4]);
hold off;
subplot(2,1,2);
plot(dt*forces(:,1), forces(:,3), plot_options)
ylabel('$C_Y$');
xlabel('$t$');
hold off;