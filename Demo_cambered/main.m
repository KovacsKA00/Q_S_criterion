clc
clear
close all

tic

c = 0.1; % chord-length of the cambered plate
res_factor = 2; % factor to change the resolution of the grid

%% Read points for cambered plate profile
fileID = fopen('coordinates_cambered_geometry.txt','r');
formatSpec = '%f %f';
size_data_cambered = [2 36];
data_cambered = fscanf(fileID,formatSpec,size_data_cambered);
cambered_x = (data_cambered(1,:))';
cambered_y = (data_cambered(2,:))';

%% Read velocity data
fid = fopen('cambered_velocity_data.txt','r');
data = textscan(fid, '%f %f %f %f %f','HeaderLines',1);
fclose(fid);

nodenumber = data{1};
X_data = data{2};
Y_data = data{3};
U_data = data{4};
V_data = data{5};

% Create uniform grid
nx = res_factor * 300;  ny = res_factor * 190;
x = linspace(-0.05,0.3,nx);
y = linspace(-0.05,0.05,ny);
[xi,yi] = meshgrid(x,y);

% Interpolation to uniform grid
u_interp = scatteredInterpolant(X_data,Y_data,U_data,'linear','none');
v_interp = scatteredInterpolant(X_data,Y_data,V_data,'linear','none');

u = u_interp(xi,yi);
v = v_interp(xi,yi);

%% Q_S-criterion
min_radius = 0.002; % for the spatial coherence metric and filtering

[Q_S, Q_S_binary] = Q_S_criterion(x,y,u,v,min_radius);

figure;
pcolor(xi(1,:)/c-1,yi(:,1)/c,Q_S);
hold on;
shading interp;
axis equal tight;
set(gca,'ydir','normal');
set(gca, 'FontSize', 12);
xlabel('x/c');
ylabel('y/c');

colorbar;
ylabel(colorbar,'Q_S/(1/s)');
xlim([-1.5 2]);
xtickformat('%.1f');
ytickformat('%.1f');

l = streamslice(x/c-1, y/c, u, v, 5);
set(l,'LineWidth',0.5)
set(l,'Color','k');

% Plot vortical regions
x_Q_S = xi(Q_S_binary);
y_Q_S = yi(Q_S_binary);
scatter(x_Q_S/c-1,y_Q_S/c,7,'filled','Markeredgecolor','r','MarkerFaceColor','r')

% Plot cambered plate profile
fill(cambered_x/c-1, cambered_y/c, 'w', 'EdgeColor', 'k', 'LineWidth', 0.5);
uistack(findobj(gca,'Type','patch'),'top');

toc



