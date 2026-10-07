function writeVTKHexahedralMesh(filename, nodes, elements, nodeRegions, materialOrientations,temp_disp)
    % This function writes a hexahedral mesh to a VTK file, including material orientation vectors.
    % Inputs:
    %   filename            - name of the VTK file
    %   nodes               - Nx3 array of node coordinates
    %   elements            - Mx8 array of hexahedral element connectivity (indices to the nodes)
    %   nodeRegions         - Mx1 array of region labels for each element (optional)
    %   materialOrientations- Mx3 array of material orientation vectors for each element (optional)
    %   temp_disp           - Nx1 array of temperature distributino for each element node (optional)
    
    % Check if nodeRegions and materialOrientations are provided, if not, assign empty arrays
    if nargin < 4
        nodeRegions = []; 
        materialOrientations = [];
        temp_disp = []; 
    elseif nargin < 5
        materialOrientations = []; 
        temp_disp = [];         
    elseif nargin < 6
        temp_disp = []; 
    end

    % Open file for writing
    fid = fopen(filename, 'w');
    
    % Write VTK file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Hexahedral mesh data with material orientation\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    
    % Write points (nodes)
    numNodes = size(nodes, 1);
    fprintf(fid, 'POINTS %d float\n', numNodes);
    for i = 1:numNodes
        fprintf(fid, '%.6f %.6f %.6f\n', nodes(i, 1), nodes(i, 2), nodes(i, 3));
    end
    
    % Write hexahedral elements
    numElements = size(elements, 1);
    fprintf(fid, 'CELLS %d %d\n', numElements, numElements * 9); % 9 = 8 nodes per element + 1 for the type identifier
    for i = 1:numElements
        fprintf(fid, '8 %d %d %d %d %d %d %d %d\n', elements(i, :) - 1); % VTK uses zero-based indexing
    end
    
    % Write cell types (12 is the VTK cell type for hexahedra)
    fprintf(fid, 'CELL_TYPES %d\n', numElements);
    for i = 1:numElements
        fprintf(fid, '12\n'); % VTK_HEXAHEDRON
    end

    % Write region labels if provided
    if ~isempty(nodeRegions)
        fprintf(fid, 'CELL_DATA %d\n', size(nodeRegions, 1));
        fprintf(fid, 'SCALARS region_labels int 1\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i = 1:size(nodeRegions, 1)
            fprintf(fid, '%d\n', nodeRegions(i));
        end
    end
    
    % Write material orientation vectors if provided
    if ~isempty(materialOrientations)
        fprintf(fid, 'VECTORS material_orientations float\n');
        for i = 1:size(materialOrientations, 1)
            fprintf(fid, '%.6f %.6f %.6f\n', materialOrientations(i, 7:9));
        end
    end

    % Write material orientation vectors if provided
    if ~isempty(temp_disp)
        fprintf(fid, 'POINT_DATA %d\n', size(temp_disp, 1));
        fprintf(fid, 'SCALARS temp_value float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i = 1:size(temp_disp, 1)
            fprintf(fid, '%0.4f\n', temp_disp(i));
        end
    end

    % Close file
    fclose(fid);
    
    disp(['VTK file with material orientation written to ', filename]);
end