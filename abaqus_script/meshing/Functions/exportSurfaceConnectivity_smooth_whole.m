function exportSurfaceConnectivity_smooth_whole(grayNodes, grayFaces, whiteNodes, whiteFaces, excel_filename)
    % Function to process and export sorted surface connectivity for both "gray" and "white" matters into an Excel file
    %
    % Inputs:
    %   grayNodes - Surface node indices for gray matter
    %   grayFaces - Surface face connectivity for gray matter
    %   mediumNodes - Surface node indices for medium layer
    %   mediumFaces - Surface face connectivity for medium layer   
    %   whiteNodes - Surface node indices for white matter
    %   whiteFaces - Surface face connectivity for white matter
    %   excel_filename - The name of the Excel file to export the results
    %
    % Outputs:
    %   Excel file with two sheets: "gray_surf_connectivity" and "white_surf_connectivity"
    
    % Perform the sortSurfaceConnectivity for the entire surface
    % FacesNum = size(grayFaces,1);
    sortedGrayWholeSurface = sortSurfaceConnectivity(grayNodes, grayFaces);
    sortedWhiteWholeSurface = sortSurfaceConnectivity(whiteNodes, whiteFaces);

    %
    % % Write the whole surface data to the Excel file
    writematrix(sortedGrayWholeSurface, excel_filename, 'Sheet', 'gray_regionAll');
    writematrix(sortedWhiteWholeSurface, excel_filename, 'Sheet', 'white_regionAll');

 
    
    % % Update the first 19 regions using the sorted whole surface connectivity
    % 
    % for Idx = 0:FacesNum-2
    %     % Find the element indices of regional faces within the original whole surface (before sorting)
    %     grayElementIndices= findRegionalIndices(grayFaces{Idx + 1}, grayFaces{FacesNum});
    %     whiteElementIndices = findRegionalIndices(whiteFaces{Idx + 1}, whiteFaces{FacesNum});
    %     mediumElementIndices = findRegionalIndices(mediumFaces{Idx + 1}, mediumFaces{FacesNum});
    % 
    %     % Extract the regional connectivity using these indices from the sorted whole surface
    %     updatedGrayFaces  = sortedGrayWholeSurface(grayElementIndices, :);
    %     updatedWhiteFaces  = sortedWhiteWholeSurface(whiteElementIndices, :);
    %     updatedMediumFaces  = sortedMediumWholeSurface(mediumElementIndices, :);
    % 
    %     % % Write updated data to Excel for the corresponding regions
    %     writematrix(updatedGrayFaces, excel_filename, 'Sheet', ['gray_region' num2str(Idx)]);
    %     writematrix(updatedWhiteFaces, excel_filename, 'Sheet', ['whiteregion' num2str(Idx)]);
    %     writematrix(updatedMediumFaces, excel_filename, 'Sheet', ['middle_region' num2str(Idx)]);
    % end
    
    function [updatedSurfaceFaces] = sortSurfaceConnectivity(surfaceNodes, surfaceFaces)
        % Step 1: Sort the surface nodes and get the indices of the sorted order
        [sortedSurfaceNodes, ~] = sort(surfaceNodes);

        updatedSurfaceFaces = surfaceFaces;  % Preallocate to maintain the same size
        for i = 1:size(surfaceFaces, 1)
            for j = 1:size(surfaceFaces, 2)
                % Update each node in the face using its position in the sorted surfaceNodes
                updatedSurfaceFaces(i, j) = find(sortedSurfaceNodes == surfaceFaces(i, j), 1);
            end
        end
    end

    function elementIndices = findRegionalIndices(regionalFaces, originalWholeSurface)

        % Identify the element indices within the original whole surface
        [~, elementIndices] = ismember(regionalFaces,originalWholeSurface,'rows');

        % Filter out invalid indices (0 values indicate unmatched faces)
        elementIndices = elementIndices(elementIndices > 0);
    end
    
    % Print success message
    fprintf('Data exported successfully to %s\n', excel_filename);
end
