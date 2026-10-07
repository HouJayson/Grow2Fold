function [selectedNodeIndices, selectedElementIndices] = createWhiteMatterBoundarySet(brainPoints, whiteMatterElements, xThreshold_low, xThreshold_high, yThreshold_low, yThreshold_high, zThreshold_low, zThreshold_high)
    % Input:
    % brainPoints - Nx3 matrix of node coordinates
    % whiteMatterElements - Mx8 matrix of element connectivity (assuming hexahedral elements)
    % xThreshold - Threshold for x-coordinate (e.g., [-0.3,0.3])
    
    % Output:
    % selectedNodeIndices - Indices of nodes where x > xThreshold
    % selectedElementIndices - Indices of elements where all nodes' x > xThreshold
    uniquePointIndices = unique(whiteMatterElements(:));


    % Find the nodes where the x-coordinate is greater than the threshold
    NodeIndices_roi = find(xThreshold_low <= brainPoints(:, 1) & brainPoints(:, 1) <= xThreshold_high & ...
                               yThreshold_low <= brainPoints(:, 2) & brainPoints(:, 2) <= yThreshold_high & ...
                               zThreshold_low <= brainPoints(:, 3) & brainPoints(:, 3) <= zThreshold_high);

    selectedNodeIndices = intersect(uniquePointIndices, NodeIndices_roi);
    
    % Initialize an empty array to store the selected elements
    selectedElementIndices = [];
    
    % Loop through each element to check if all nodes satisfy the condition
    for i = 1:size(whiteMatterElements, 1)
        elementNodes = whiteMatterElements(i, :);
        if all(ismember(elementNodes, selectedNodeIndices))
            selectedElementIndices = [selectedElementIndices; i];
        end
    end
end
