% =============================================================================
% Copyright (c) 2026 Jixin Hou et al.
% All rights reserved.
%
% This code is provided as part of the research software accompanying:
% [Grow2Fold: mapping heterogeneous developmental growth to human brain folding]
%
% Use, modification, and redistribution are permitted under the terms of the
% license provided in the LICENSE file of this repository.
%
% Repository: [https://github.com/BioDMX-UGA/Grow2Fold]
% =============================================================================

clc; clear;

addpath('E:\PhD\PINN\WholeBrain_Pipeline\postprocessing\functionsMat');

currentDir   = "E:\PhD\PINN\WholeBrain_Pipeline\postprocessing";
fPath_ref    = fullfile(currentDir, 'datafromAbaqus','refBrain');
fPath_conn   = fullfile(currentDir, 'datafromAbaqus','connFile');
fPath_input  = fullfile(currentDir, 'datafromAbaqus','rawData');

caseId       = {'CC01022XX08'};
conditions   = 'par_Clustering';
casename     = 'Hex_thk14_skull_2steps_stress';

for idx = 1:length(caseId)
    dataStoragePath    = fullfile(currentDir,'dataSummary',caseId{idx},'whole');
    if ~exist(dataStoragePath, 'dir')
        mkdir(dataStoragePath);
    end
    num_frames   = 26;   
    frame_interval = 1;
    iterations   = 10;
    lambda       = 0.2;
    batchSize    = 5000;

    sheets = struct('gray','gray_regionAll','white','white_regionAll');

    % === Read invariant files once ===
    connectivity_file    = fullfile(fPath_conn, sprintf('Sconn_%s.xlsx', caseId_ref{idx}));
    reference_file_inner = fullfile(fPath_ref, sprintf('%s.InnerSurf_W.vtk', caseId_ref{idx}));
    reference_file_outer = fullfile(fPath_ref, sprintf('%s.OuterSurf_W.vtk', caseId_ref{idx}));

    SurfConnectivity.gray = readmatrix(connectivity_file, 'Sheet', sheets.gray);
    SurfConnectivity.white = readmatrix(connectivity_file, 'Sheet', sheets.white);

    pvtk_inner = mvtk_read(reference_file_inner);
    pvtk_outer = mvtk_read(reference_file_outer);

    par_info_huang  = pvtk_inner.par_huang;
    par_info_FS2009 = pvtk_inner.par_FS2009;
    par_info_Clustering  = pvtk_inner.par_Clustering;

    Nodes_par.inner = pvtk_inner.vertices;
    Nodes_par.outer = pvtk_outer.vertices;

    kdtree.inner = createns(Nodes_par.inner, 'NSMethod', 'kdtree');
    kdtree.outer = createns(Nodes_par.outer, 'NSMethod', 'kdtree');

    fprintf('Processing case: %s\n', caseId{idx});

    coordinate_file_gray = fullfile(fPath_input, sprintf('%s.gray.%s.%s.mat', caseId{idx}, conditions, casename) );
    coordinate_file_white = fullfile(fPath_input, sprintf('%s.white.%s.%s.mat', caseId{idx}, conditions, casename) );

    S_gray = load(coordinate_file_gray);
    S_white = load(coordinate_file_white);
    grayData = S_gray.data;
    whiteData = S_white.data;
    NodeCoord_gray0 = grayData(:,1:3,1);
    NodeCoord_white0 = whiteData(:,1:3,1);

    triConn_gray = quadToTriConnectivity(SurfConnectivity.gray, NodeCoord_gray0);
    triConn_white = quadToTriConnectivity(SurfConnectivity.white, NodeCoord_white0);

    smoothed_gray0 = laplacianSmoothing(triConn_gray, NodeCoord_gray0, iterations, lambda);
    smoothed_white0 = laplacianSmoothing(triConn_white, NodeCoord_white0, iterations, lambda);

    % === Assign parcellation ===
    par_simulation_huang_gray  = assignParcellation(smoothed_gray0, kdtree.outer, par_info_huang, batchSize);
    par_simulation_FS2009_gray = assignParcellation(smoothed_gray0, kdtree.outer, par_info_FS2009, batchSize);
    par_simulation_Clustering_gray = assignParcellation(smoothed_gray0, kdtree.outer, par_info_Clustering, batchSize);

    par_simulation_huang_white  = assignParcellation(smoothed_white0, kdtree.inner, par_info_huang, batchSize);
    par_simulation_FS2009_white = assignParcellation(smoothed_white0, kdtree.inner, par_info_FS2009, batchSize);
    par_simulation_Clustering_white = assignParcellation(smoothed_white0, kdtree.inner, par_info_Clustering, batchSize);   

    % === Process frames in parallel ===
    parfor i = 1:num_frames
        NodeCoord_gray = grayData(:,1:3,i);
        NodeCoord_white = whiteData(:,1:3,i);
        smoothed_gray = laplacianSmoothing(triConn_gray, NodeCoord_gray, iterations, lambda);
        smoothed_white = laplacianSmoothing(triConn_white, NodeCoord_white, iterations, lambda);
        vtkFilename_gray = char(fullfile(dataStoragePath, sprintf('%s_gray_frame%d.vtk', caseId{idx}, frame_interval*(i-1))));
        vtkFilename_white = char(fullfile(dataStoragePath, sprintf('%s_white_frame%d.vtk', caseId{idx}, frame_interval*(i-1))));

        dataVTK_gray = struct('vertices',smoothed_gray,'faces',triConn_gray, ...
            'par_huang',par_simulation_huang_gray, 'par_FS2009',par_simulation_FS2009_gray,'par_Clustering', par_simulation_Clustering_gray);
        dataVTK_white = struct('vertices',smoothed_white,'faces',triConn_white, ...
            'par_huang',par_simulation_huang_white, 'par_FS2009',par_simulation_FS2009_white,'par_Clustering', par_simulation_Clustering_white);

        mvtk_write(dataVTK_gray, vtkFilename_gray, 'legacy');
        mvtk_write(dataVTK_white, vtkFilename_white, 'legacy');
    end

    GI = zeros(num_frames,1);
    GI_alpha = zeros(num_frames,1);

    for i=1:num_frames  
        nvtk_gray = fullfile(dataStoragePath, sprintf('%s_gray_frame%d.vtk',caseId{idx}, i-1));
        pvtk_gray=mvtk_read(nvtk_gray);
        nvtk_white = fullfile(dataStoragePath, sprintf('%s_white_frame%d.vtk',caseId{idx}, i-1));
        pvtk_white=mvtk_read(nvtk_white);  

        Nodes_gray=pvtk_gray.vertices;
        Nodes_white=pvtk_white.vertices;

        connectivity_gray=double(pvtk_gray.faces);

        Cthickness = corticalThickness(Nodes_gray, Nodes_white);

        Cthickness1 = dataSmoothing(Nodes_gray, connectivity_gray, Cthickness, 2, 1);

        Cthickness1(pvtk_white.par_huang==0)=0;

        [GC, MC_dimensionless, MC, k_max, k_min, shapeIndex]= Curvatures(connectivity_gray,Nodes_gray(:,1),Nodes_gray(:,2),Nodes_gray(:,3));  
        data_to_smooth = [GC, MC_dimensionless, MC, k_max, k_min, shapeIndex];
        smoothed_data = dataSmoothing(Nodes_gray, connectivity_gray, data_to_smooth, 10, 0.5);

        GC_corr = smoothed_data(:,1);
        MC_corr = smoothed_data(:,2);
        MC_dimensionless_corr = smoothed_data(:,3);
        k_max_corr = smoothed_data(:,4);
        k_min_corr = smoothed_data(:,5);
        shapeIndex_corr = smoothed_data(:,6);

        [Connectivity_hull, vertices_hull, GI(i)]= GyrificationIndex(connectivity_gray,Nodes_gray(:,1),Nodes_gray(:,2),Nodes_gray(:,3)); 
        SulcDepth = SulcalDepth(Connectivity_hull,vertices_hull,Nodes_gray(:,1),Nodes_gray(:,2),Nodes_gray(:,3));  

        [Connectivity_hull_alpha, vertices_hull_alpha, GI_alpha(i)]= GyrificationIndex_alphashape(connectivity_gray,Nodes_gray(:,1),Nodes_gray(:,2),Nodes_gray(:,3));  
        SulcDepth_alpha = SulcalDepth_alphashape(Connectivity_hull_alpha,vertices_hull_alpha,Nodes_gray(:,1),Nodes_gray(:,2),Nodes_gray(:,3));  
        SulcDepth_alpha_shrink = SulcalDepth_alphashape_shrink(Connectivity_hull_alpha,vertices_hull_alpha,Nodes_gray(:,1),Nodes_gray(:,2),Nodes_gray(:,3));  

        data_to_smooth1= [SulcDepth,SulcDepth_alpha,SulcDepth_alpha_shrink];
        smoothed_data1 = dataSmoothing(Nodes_gray, connectivity_gray, data_to_smooth1, 2, 1);

        SulcDepth_corr = smoothed_data1(:,1);
        SulcDepth_alpha_corr = smoothed_data1(:,2);
        SulcDepth_alpha_shrink_corr = smoothed_data1(:,3);

        dataVTK = pvtk_gray;
        dataVTK.MC=MC_corr;
        dataVTK.Cthickness=Cthickness1;
        dataVTK.MC_dimensionless=MC_dimensionless_corr;
        dataVTK.k_max=k_max_corr;
        dataVTK.k_min=k_min_corr;
        dataVTK.shapeIndex=shapeIndex_corr;
        dataVTK.SulcDepth=SulcDepth_corr;
        dataVTK.SulcDepth_alpha=SulcDepth_alpha_corr;
        dataVTK.SulcDepth_alpha_shrink=SulcDepth_alpha_shrink_corr;
        dataVTK.GI=GI(i)*ones(length(MC_corr),1);
        dataVTK.GI_alpha=GI_alpha(i)*ones(length(MC_corr),1);

        mvtk_write(dataVTK,char(nvtk_gray),'legacy');
        fprintf('Frame %d, finished!\n',(i-1));
    end
end
disp('All vtk files have been written.');

function par_cell = assignParcellation(nodes, kdtree, par_info, batchSize)
    numPoints  = size(nodes, 1);
    numBatches = ceil(numPoints / batchSize);
    par_cell   = cell(numBatches, 1);

    for batch = 1:numBatches
        idxRange = (batch - 1) * batchSize + 1 : min(batch * batchSize, numPoints);
        batchPoints = nodes(idxRange, :);
        idx = knnsearch(kdtree, batchPoints);
        par_cell{batch} = par_info(idx);
    end

    par_cell = vertcat(par_cell{:});
end