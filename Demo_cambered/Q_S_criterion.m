function [Q_S, Q_S_binary] = Q_S_criterion(x,y,u,v,min_radius_phys)
% Compute gradients
[ux,uy] = gradient(u,x,y);
[vx,vy] = gradient(v,x,y);

% Compute Q-criterion
Q = zeros(size(u));
for i = 1:size(u,1)
    for j = 1:size(u,2)
        G = [ux(i,j), uy(i,j);
             vx(i,j), vy(i,j)];
        Q_S = 0.5*(G + G');     % strain rate tensor
        Omega = 0.5*(G - G');   % rotation tensor
        Q(i,j) = 0.5*(norm(Omega,'fro')^2 - norm(Q_S,'fro')^2);
    end
end

% Vorticity
omega = vx - uy;
omega = imgaussfilt(omega,1); % smoothing

D = nan(size(omega));

% Compute radius in grid units
[Ny, Nx] = size(omega);
dx = abs(x(2) - x(1));
dy = abs(y(2) - y(1));
r_x = ceil(min_radius_phys / dx);
r_y = ceil(min_radius_phys / dy);

% Create circular mask based on physical distance
[X, Y] = meshgrid((-r_x:r_x)*dx, (-r_y:r_y)*dy);
mask = (X.^2 + Y.^2) <= min_radius_phys^2;

% Pad omega to avoid border issues
omega_padded = padarray(omega, [r_y, r_x], 'replicate');

for i = 1:Ny
    for j = 1:Nx
        % Define local window
        win = omega_padded(i:i+2*r_y, j:j+2*r_x);
        % Apply mask to extract values
        values = win(mask);
        D(i,j) = std(values);
    end
end

% Coherence weight
C = 1./D;

% Combine Q and coherence
Q_S = Q.*C;

% Vortex identification using GMM thresholding
Q_S_nonzero = Q_S(Q_S > 0); % select positive Q_S values
options = statset('MaxIter',1000);
GMModel = fitgmdist(Q_S_nonzero, 2, 'Options', options); % two Gaussian distributions to model the background, and the vortical regions

mu1 = GMModel.mu(1); mu2 = GMModel.mu(2); % mean
sigma1 = sqrt(GMModel.Sigma(:,:,1)); sigma2 = sqrt(GMModel.Sigma(:,:,2)); % variance
p1 = GMModel.ComponentProportion(1); p2 = GMModel.ComponentProportion(2); % prior probabilities

pdf1 = @(x) p1 * normpdf(x, mu1, sigma1);
pdf2 = @(x) p2 * normpdf(x, mu2, sigma2);
f = @(x) pdf1(x) - pdf2(x);
Q_S_iso = fzero(f, mean([mu1 mu2])); % finds the value of x where the two Gaussians are equally likely

% Create binary mask and filter by area
Q_S_binary = Q_S > Q_S_iso;
grid_point_area = abs((x(2)-x(1)) * (y(2)-y(1))); % in m^2
min_area_phys = min_radius_phys^2*pi; % minimum area in m^2
min_area_grid_points = round(min_area_phys / grid_point_area);  % minimum number of grid points in a connected region

% Filter small regions
% Remove all connected components from Q_S_binary that have fewer grid points than min_area_grid_points
Q_S_binary = bwareaopen(Q_S_binary, min_area_grid_points);


end
