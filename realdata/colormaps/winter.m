function map = winter(siz)
% Initialize RGB nodes
rgb_nodes = [
    0   0   255;
    0   128 196;
    0   255 128];

% Normalize to [0, 1]
rgb_nodes = rgb_nodes / 255.0;

% Initialize output colormap
map = zeros(siz, 3);

% Interpolate RGB
for i = 1:3
    map(:,i) = interp1(linspace(0,1,3),rgb_nodes(:,i),linspace(0,1,siz));
end
