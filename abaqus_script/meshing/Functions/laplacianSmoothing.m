function smoothedNodes = laplacianSmoothing(surfaceFaces, brainNodes, iterations, lambda)
    % surfaceFaces: Connectivity of the surface (face indices)
    % brainNodes: Coordinates of the mesh nodes
    % iterations: Number of smoothing iterations
    % lambda: Smoothing factor (usually between 0 and 1)

    % Initialize the smoothed nodes as the original nodes
    smoothedNodes = brainNodes;

    % Precompute the neighbors for each node based on surface connectivity
    neighbors = cell(size(brainNodes, 1), 1); % Create a cell array to hold neighbors for each node
    for i = 1:size(surfaceFaces, 1)
        for j = 1:size(surfaceFaces, 2)
            current_node = surfaceFaces(i, j);
            neighbors{current_node} = unique([neighbors{current_node}; surfaceFaces(i, :)']);
        end
    end

    % Remove self from neighbors
    for i = 1:length(neighbors)
        if ~isempty(neighbors{i})
            neighbors{i}(neighbors{i} == i) = []; % Remove self if present
        end
    end

    % Perform the smoothing iterations
    for iter = 1:iterations
        % Compute new coordinates for all nodes simultaneously
        newCoords = smoothedNodes; % Initialize with current coordinates

        % Update node positions based on neighbors
        parfor i = 1:size(brainNodes, 1)
            if ~isempty(neighbors{i})
                % Calculate the centroid of neighboring nodes
                centroid = mean(smoothedNodes(neighbors{i}, :), 1);

                % Update the position of the current node (weighted towards centroid)
                newCoords(i, :) = smoothedNodes(i, :) + lambda * (centroid - smoothedNodes(i, :));
            end
        end

        % Update the node positions after all updates are computed
        smoothedNodes = newCoords;
    end
end
