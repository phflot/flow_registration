% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.
    
function quiver_viz = quiver_visualization(img, w, scale)
%QUIVER_VISULAIZATION function for the quiver visualization
%   Takes an image and a displacementfield w as input

    if nargin < 3
        scale = 1;
    end

    [n, m, ~] = size(w);

    w = imresize(w, 0.03);    
    X_w = linspace(1, m, size(w, 2));
    Y_w = linspace(1, n, size(w, 1));
    [X_m, Y_m] = meshgrid(X_w, Y_w);
    X_m = X_m + 0.5 * (X_w(2) - X_w(1));
    Y_m = Y_m + 0.5 * (Y_w(2) - Y_w(1));
    
    imshow(img);
    hold on
    lh = streamline(X_w, Y_w, w(:, :, 1), w(:, :, 2), X_m, Y_m);
    
    X_w = X_w(:, 2:end-1);
    Y_w = Y_w(:, 2:end-1);
    w = w(2:end-1, 2:end-1, :);
    quiver(X_w, Y_w, w(:, :, 1), w(:, :, 2), scale, 'LineWidth', 2, ...
    'Color', 'w');
    for i = 1:length(lh)
        lh(i).Color=[0,0,0,0.75];
    end
    quiver_viz = getframe(gca).cdata;
    set(gca,'visible','off')
    hold off
end

