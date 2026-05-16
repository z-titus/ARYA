%%%%%%%%%%%%
% build quiver plots of output SIMPLE data
%%%%%%%%%%
save_vector = 0; % set to 1 to save plot
save_stream = 0;
save_p = 0;
save_res = 0;
Nx = 100;   
Ny = 100;   

Lx = 1.0; 
Ly = 1.0;  

path_prfx = '<PATH/TO/ARYA>';
path = '<PATH/TO/OUTPUT/DATA>';
file_u = strcat(path_prfx,  path,  '\data_u');
file_v = strcat(path_prfx, path, '\data_v');
file_p = strcat(path_prfx, path, '\data_p');
file_res = strcat(path_prfx, path, '\data_res');


fid = fopen(file_u,'r');
u_data = fscanf(fid,'%f');
fclose(fid);

fid = fopen(file_v,'r');
v_data = fscanf(fid,'%f');
fclose(fid);

fid = fopen(file_p,'r');
p_data = fscanf(fid,'%f');
fclose(fid);
p = reshape(p_data, [Nx+2, Ny+2])';

u = reshape(u_data, [Nx+1, Ny+2])';
v = reshape(v_data, [Nx+2, Ny+1])';



[u_c, v_c] = center_data(u, v, Nx, Ny);

% x = [0, linspace(Lx/(2*Nx), Lx - Lx/(2*Nx), Nx), Lx];
% y = [0, linspace(Ly/(2*Ny), Ly - Ly/(2*Ny), Ny), Ly];

x = linspace(Lx/(2*Nx), Lx - Lx/(2*Nx), Nx);
y = linspace(Ly/(2*Ny), Ly - Ly/(2*Ny), Ny);

[X, Y] = meshgrid(x, y);

% Uplot = [zeros(Ny+2, 1), u_c, zeros(Ny+2, 1)];
% Vplot = [zeros(1, Nx+2); v_c; zeros(1, Nx+2)];

Uplot = u_c(2:end-1,:);
Vplot = v_c(:,2:end-1);

plot_dens = 1;

f = figure(1001);
quiver(X(1:plot_dens:end, 1:plot_dens:end), ...
       Y(1:plot_dens:end, 1:plot_dens:end), ...
       Uplot(1:plot_dens:end, 1:plot_dens:end), ...
       Vplot(1:plot_dens:end, 1:plot_dens:end), 5);
hold on
xlabel('$x_1$', 'Fontsize', 24, 'Interpreter', 'latex');
ylabel('$x_2$','Fontsize', 24, 'Interpreter', 'latex');
xlim([-0.05 1.05])
ylim([-0.05 1.05])
%title('Vector Field - $Re = 100$','Fontsize', 16, 'Interpreter', 'latex');
axis equal tight;
if save_vector == 1
    %plot(0.566, 0.613, marker='.', color='red', markersize=16)
    %set(gcf, 'Position', [100 100 600 600])
    theme(f, "light")
    saveas(gcf,'vctr_fld_Re400_100-100','epsc')
end 
hold off

stream_dens = 1;
[xs1, ys1] = meshgrid(linspace(0.0,0.1,3), linspace(0.0,0.07,3));
[xs2, ys2] = meshgrid(linspace(0.7,0.95,3), linspace(0.05,0.3,3));

f = figure(2001);
streamline(X, Y, Uplot, Vplot, xs1, ys1, 'Color', [0 0.4470 0.7410])
hold on
streamline(X, Y, Uplot, Vplot, xs2, ys2, 'Color', [0 0.4470 0.7410])
hold on
streamline(X, Y, Uplot, Vplot, X(1:95:end), Y(1:95:end), 'Color', [0 0.4470 0.7410])
hold off

xlabel('$x_1$', 'Fontsize', 24, 'Interpreter', 'latex');
ylabel('$x_2$','Fontsize', 24, 'Interpreter', 'latex');
xlim([-0.05 1.05])
ylim([-0.05 1.05])
%title('Streamlines','Fontsize', 16, 'Interpreter', 'latex');
axis equal tight;
if save_stream == 1
    theme(f, "light")
    set(gcf, 'Position', [100 100 600 600]) % pixels
    exportgraphics(gcf, 'streamlines_Re400_100-100.eps')
end 

p_plt = p(2:Nx+1, 2:Ny+1);

f = figure(3001);
contourf(X,Y, p_plt, 100, 'linecolor', 'None')
daspect([1 1 1])
c = colorbar();
colormap('turbo')
cL = ylabel(c, 'p''', 'fontsize',24, 'interpreter','latex');
xlabel('$x_1$', 'Fontsize', 24, 'Interpreter', 'latex');
ylabel('$x_2$','Fontsize', 24, 'Interpreter', 'latex');
set(cL, 'Rotation', 90)
if save_p == 1
    theme(f, "light")
    %set(gcf, 'Position', [100 100 600 600]) % pixels
    exportgraphics(gcf, 'pressure_Re100_100-100.png', 'Resolution',300)
end 

res_data = readmatrix(file_res, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);

if length(size(res_data))>1
    res_data = res_data(:,2);
end 

f = figure(4001);
semilogy(linspace(5,length(res_data),length(res_data)-4), res_data(5:end)/res_data(5), marker= "+");
xlabel('Time step', 'Fontsize', 20, 'Interpreter', 'latex');
ylabel('$\frac{\Sigma b''_{ij}}{\Sigma b''_{ij}(5)}$', 'Fontsize', 20, 'Interpreter', 'latex');
if save_res == 1
    theme(f, "light")
    exportgraphics(gcf, 'res_Re1000_100-100.eps')
end 


function [u_c, v_c] = center_data(u, v, Nx, Ny)

    % average u in x-direction
    u_c = 0.5 * (u(:, 1:Nx) + u(:, 2:Nx+1));

    % average v in y-direction
    v_c = 0.5 * (v(1:Ny, :) + v(2:Ny+1, :));

end
% validation
figure(5001)
x1 = [0.614 0.555 0.531];
x2 = [0.732 0.614 0.571];

a = plot(x1,x2, marker = 'o');
labels = ['Re = 100', 'Re = 400', 'Re = 1000'];

%% validation
save_val = 0;
f = figure(6001);
x1 = [0.6172 0.5547 0.5313];
x2 = [0.7344 0.6055 0.5625];

b = plot(x1,x2, marker = 'o', MarkerSize=10);
text(x1-0.017, x2+0.015, {'$Re = 100$','$Re = 400$','$Re = 1000$'}, 'interpreter', 'latex', 'FontSize', 13)
hold on
% x1 = [0.62, 0.57];
% x2 = [0.74, 0.62];
x1 = [0.621, 0.566, 0.545];
x2 = [0.738, 0.613, 0.565];
xlim([0.5 0.64])
ylim([0.4 0.9])
plot(x1,x2, marker = '+', MarkerSize=12);
xlabel('$x_1$', 'Fontsize', 24, 'Interpreter', 'latex');
ylabel('$x_2$','Fontsize', 24, 'Interpreter', 'latex');
lgd = legend('Ghia et al.', 'Current study', 'Interpreter', 'latex', 'Fontsize', 13, 'Location','northwest');
lgd.Position(3) = lgd.Position(3) * 1.2;  % widen legend box
hold off
if save_val == 1
    theme(f, "light")
    %set(gcf, 'Position', [100 100 600 600]) % pixels
    exportgraphics(gcf, 'validation.eps')
end 
%% internal geometry
save_vector=0;
save_stream=0;
save_h=0;

Nx = 100;   
Ny = 100;   

% internal box geometry
Nx_start=40;
Nx_end=60;
Ny_start=40;
Ny_end=60;

Lx = 1.0; 
Ly = 1.0;  

path_prfx = '\\wsl.localhost\Ubuntu\home\zach_t\me6313\ARYA\_case';
path = '\sliding_lid\int_geometry';
file_u = strcat(path_prfx,  path,  '\data_u');
file_v = strcat(path_prfx, path, '\data_v');

fid = fopen(file_u,'r');
u_data = fscanf(fid,'%f');
fclose(fid);

fid = fopen(file_v,'r');
v_data = fscanf(fid,'%f');
fclose(fid);

u = reshape(u_data, [Nx+1, Ny+2])';
v = reshape(v_data, [Nx+2, Ny+1])';

[u_c, v_c] = center_data(u, v, Nx, Ny);

x = linspace(Lx/(2*Nx), Lx - Lx/(2*Nx), Nx);
y = linspace(Ly/(2*Ny), Ly - Ly/(2*Ny), Ny);

[X, Y] = meshgrid(x, y);

Uplot = u_c(2:end-1,:);
Vplot = v_c(:,2:end-1);

plot_dens = 2;

f = figure(1001);
quiver(X(1:plot_dens:end, 1:plot_dens:end), ...
       Y(1:plot_dens:end, 1:plot_dens:end), ...
       Uplot(1:plot_dens:end, 1:plot_dens:end), ...
       Vplot(1:plot_dens:end, 1:plot_dens:end), 6.5);
hold on
xlabel('$x_1$', 'Fontsize', 24, 'Interpreter', 'latex');
ylabel('$x_2$','Fontsize', 24, 'Interpreter', 'latex');
xlim([-0.05 1.05])
ylim([-0.05 1.05])
%title('Vector Field - $Re = 100$','Fontsize', 16, 'Interpreter', 'latex');
mask = false(Ny, Nx);
mask(Ny_start:Ny_end, Nx_start:Nx_end) = true;
contourf(X-Lx/(2*Nx), Y+Ly/(2*Ny), mask, [1 1], 'FaceColor', 'red', 'FaceAlpha', 0.7, 'LineStyle', 'none');
hold off
axis equal tight;
if save_vector == 1
    theme(f, "light")
    set(gcf, 'Position', [100 100 600 600]) % pixels
    saveas(gcf,'vctr_fld_Re400_int_g','epsc')
end 

stream_dens = 75;
[xs1, ys1] = meshgrid(linspace(0.0,0.1,3), linspace(0.0,0.07,3));
[xs2, ys2] = meshgrid(linspace(0.7,0.95,3), linspace(0.05,0.3,3));
[xs3,ys3]  = meshgrid(linspace(0.35,0.6,8), linspace(0.35,0.65,8));

f = figure(2001);
streamline(X, Y, Uplot, Vplot, xs1, ys1, 'Color', [0 0.4470 0.7410])
hold on
streamline(X, Y, Uplot, Vplot, xs2, ys2, 'Color', [0 0.4470 0.7410])
hold on
streamline(X, Y, Uplot, Vplot, xs3, ys3, 'Color', [0 0.4470 0.7410])
hold on
streamline(X, Y, Uplot, Vplot, X(1:stream_dens:end), Y(1:stream_dens:end), 'Color', [0 0.4470 0.7410])
mask = false(Ny, Nx);
mask(Ny_start:Ny_end, Nx_start:Nx_end) = true;
contourf(X-Lx/(Nx), Y+Ly/(2*Ny), mask, [1 1], 'FaceColor', 'red', 'FaceAlpha', 0.7, 'LineStyle', 'none');
hold off
xlabel('$x_1$', 'Fontsize', 24, 'Interpreter', 'latex');
ylabel('$x_2$','Fontsize', 24, 'Interpreter', 'latex');
xlim([-0.05 1.05])
ylim([-0.05 1.05])
%title('Streamlines','Fontsize', 16, 'Interpreter', 'latex');
axis equal tight;
if save_stream == 1
    theme(f, "light")
    set(gcf, 'Position', [100 100 600 600]) % pixels
    exportgraphics(gcf, 'streamlines_Re400_int_g.eps')
end 



file_h = strcat(path_prfx,  path,  '\data_h');

fid = fopen(file_h,'r');
h_data = fscanf(fid,'%f');
fclose(fid);

h = reshape(h_data, [Nx+2, Ny+2]);

h_plt = h(2:Nx+1, 2:Ny+1)';

x = linspace(Lx/(2*Nx), Lx - Lx/(2*Nx), Nx);
y = linspace(Ly/(2*Ny), Ly - Ly/(2*Ny), Ny);

[X, Y] = meshgrid(x, y);

f = figure(7001);
contourf(X,Y, h_plt, 100, 'linecolor', 'None')
daspect([1 1 1])
c = colorbar();
c.FontSize = 14;
%clim([-1 50]); 
colormap('turbo')
cL = ylabel(c, 'h (kJ/kg)', 'fontsize',20, 'interpreter','latex');
set(cL, 'Rotation', 90)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14; 
xlabel('$x_1$', 'Fontsize', 24, 'Interpreter', 'latex');
ylabel('$x_2$','Fontsize', 24, 'Interpreter', 'latex');
if save_h == 1
    theme(f, "light")
    exportgraphics(gcf, 'energy_Re400_int_g_t180.png', 'Resolution',300)              
end 

