function out = build_theta_indicator_temperature( ...
        brainPoints, whiteElements, ...
        pvertice_gray, pface_gray, pvertice_white, ...
        pvertice_gray_mesh, ...
        par_MMP, par_FS2009, ...
        iterations, alphaRadius, cortical_thickness)
%BUILD_THETA_INDICATOR_TEMPERATURE
% Build indicator masks and apply temperature distribution theta_y.
%
% Inputs
%   brainPoints          (nNodes x 3) FE mesh node coordinates
%   whiteElements        (nWhiteElem x nen) connectivity (node IDs)
%   pvertice_gray        (nV x 3) pial/gray surface vertices
%   pface_gray           (nF x 3) pial/gray surface faces (1-based)
%   pvertice_white       (nV x 3) white surface vertices (same indexing as gray)
%   pvertice_gray_mesh   (nVm x 3) vertices used to build KDTree for distance
%   par_MMP              (nV x 1) MMP labels per surface vertex
%   par_FS2009           (nV x 1) FS2009 labels per surface vertex
%   iterations           (1 x nROI) shrink iterations for each ROI (same order as ROI list below)
%   alphaRadius          alphaShape radius (e.g., 0.9)
%   cortical_thickness   scalar (e.g., 1.0)
%
% Output struct out:
%   out.theta_y          (nNodes x 1)
%   out.theta_y_temp     (nNodes x 1)
%   out.insideMask_sum   (nNodes x 1) logical
%   out.distToGraySurface(nNodes x 1)
%
% Requires: inpolyhedron, alphaShape toolbox

    if nargin < 9 || isempty(iterations)
        error('iterations must be provided (1 x 14 in your current setup).');
    end
    if nargin < 10 || isempty(alphaRadius), alphaRadius = 1.5; end
    if nargin < 11 || isempty(cortical_thickness), cortical_thickness = 1.4; end

    % Ensure faces are 1-based
    if min(pface_gray(:)) == 0
        pface_gray = pface_gray + 1;
    end

    % -------- 1) Distance-based baseline theta_y_temp --------
    kdtree = createns(pvertice_gray_mesh, 'NSMethod', 'kdtree');
    [~, distToGraySurface] = knnsearch(kdtree, brainPoints);

    theta_y_temp = (1 ./ (1 + exp(10 .* (distToGraySurface ./ cortical_thickness - 1))) + 0.1) ./ 1.1;
    theta_y = theta_y_temp;

    % -------- 2) Define ROIs (exactly as your code) --------
    % roi_indices1  = ismember(par_MMP,[1356321,1428771,1497135]);                 % central sulcus
    % roi_indices2  = ismember(par_MMP,2539284);                                   % postcentral sulcus
    % roi_indices3  = ismember(par_FS2009,3988703);                                % superior temporal sulcus
    % roi_indices4  = ismember(par_MMP,2137626);                                   % precentral sulcus
    % roi_indices5  = ismember(par_MMP,[8691849,6256510,5726817]);                 % inferior frontal sulcus
    % roi_indices6  = ismember(par_FS2009,6609981);                                % superior frontal sulcus
    % roi_indices7  = ismember(par_MMP,[10935432,10383679,13284507, ...
    %                                  13422216,13607024,11445936,8419702]);      % intraparietal sulcus
    % 
    % roi_indices8  = ismember(par_MMP,[1751584,1958189,1623337]);                 % (right) central sulcus
    % roi_indices9  = ismember(par_MMP,2934557);                                   % (right) postcentral sulcus
    % roi_indices10 = ismember(par_FS2009,3988704);                                % (right) superior temporal sulcus
    % roi_indices11 = ismember(par_MMP,2729248);                                   % (right) precentral sulcus
    % roi_indices12 = ismember(par_MMP,[8952963,6847371,7303293]);                 % (right) inferior frontal sulcus
    % roi_indices13 = ismember(par_FS2009,6609982);                                % (right) superior frontal sulcus
    % roi_indices14 = ismember(par_MMP,[11330693,11908538,13089682, ...
    %                                  12964479,8485757]);                        % (right) intraparietal sulcus

    roi_indices1  = ismember(par_FS2009,660701);                        % (left) central sulcus
    roi_indices2  = ismember(par_FS2009,[9221140,13143061]);            % (left) postcentral sulcus
    roi_indices3  = ismember(par_FS2009,3988703);                       % (left) superior temporal sulcus
    roi_indices4  = ismember(par_FS2009,[13112341,15733781]);           % (left) precentral sulcus
    roi_indices5  = ismember(par_FS2009,1367261);                       % (left) inferior frontal sulcus
    roi_indices6  = ismember(par_FS2009,[6609981,6558861]);             % (left) superior frontal sulcus
    roi_indices7  = ismember(par_FS2009,14423183);                      % (left) intraparietal sulcus
    % roi_indices8  = ismember(par_FS2009,11842623);                      % (left) parieto-occipital & Calcarine sulcus
    % roi_indices9  = ismember(par_FS2009,[1346781,9180220,13158500]);    % (left) collateral sulcus
    % roi_indices10  = ismember(par_FS2009,[6558941,9845786,...
                                   % 4930586,15386]);                     % (left) cingulate sulcus

    roi_indices11  = ismember(par_FS2009,660702);                       % (right) central sulcus
    roi_indices12  = ismember(par_FS2009,[9221141,13143062]);           % (right) postcentral sulcus
    roi_indices13 = ismember(par_FS2009,3988704);                       % (right) superior temporal sulcus
    roi_indices14 = ismember(par_FS2009,[13112342,15733782]);           % (right) precentral sulcus
    roi_indices15 = ismember(par_FS2009,1367262);                       % (right) inferior frontal sulcus
    roi_indices16 = ismember(par_FS2009,[6609982,6558862]);             % (right) superior frontal sulcus
    roi_indices17 = ismember(par_FS2009,14423184);                      % (right) intraparietal sulcus
    % roi_indices18  = ismember(par_FS2009,11842624);                     % (right) parieto-occipital & Calcarine sulcus
    % roi_indices19  = ismember(par_FS2009,[1346782,9180221,13158501]);   % (right) collateral sulcus
    % roi_indices20  = ismember(par_FS2009,[6558942,9845787,...
                                   % 4930587,15387]);                     % (right) cingulate sulcus

    roi_indices_sum = [roi_indices1,roi_indices2,roi_indices3,roi_indices4,roi_indices5,roi_indices6,roi_indices7,...
                       roi_indices11,roi_indices12,roi_indices13,roi_indices14,roi_indices15,roi_indices16,roi_indices17];

    nROI = size(roi_indices_sum,2);
    if numel(iterations) ~= nROI
        error('iterations length (%d) must match number of ROIs (%d).', numel(iterations), nROI);
    end

    % -------- 3) Build insideMask_sum by alphaShape + inpolyhedron --------
    insideMask_sum = false(size(brainPoints,1), 1);

    for j = 1:nROI
        roi_indices = roi_indices_sum(:,j);
        roi_vertices = find(roi_indices);

        if isempty(roi_vertices)
            continue;
        end

        num_iterations = iterations(j);

        % --- shrink loop (boundary stripping) ---
        for iter = 1:num_iterations
            face_in_region = all(ismember(pface_gray, roi_vertices), 2);
            region_faces1  = pface_gray(face_in_region, :);
            if isempty(region_faces1), break; end

            node_count = histcounts(region_faces1(:), 1:(size(pvertice_gray,1)+1));
            boundary_nodes = node_count > 0 & node_count < 4;

            roi_indices(boundary_nodes) = 0;
            roi_vertices = find(roi_indices);
            if isempty(roi_vertices), break; end
        end

        % --- final ROI vertices based on faces fully inside region ---
        roi_indices_final = false(size(roi_indices));
        face_in_final_region = all(ismember(pface_gray, roi_vertices), 2);
        final_region_faces = pface_gray(face_in_final_region, :);
        if isempty(final_region_faces), continue; end

        roi_indices_final(unique(final_region_faces)) = true;

        grayNodes_roi_final  = pvertice_gray(roi_indices_final, :);
        whiteNodes_roi_final = pvertice_white(roi_indices_final, :);
        if isempty(grayNodes_roi_final), continue; end

        % regionPoints (same as your current choice)
        interpolatedNodes1 = 0.5 * (grayNodes_roi_final + whiteNodes_roi_final);
        interpolatedNodes2 = 1.5 * (grayNodes_roi_final - whiteNodes_roi_final) + grayNodes_roi_final;
        interpolatedNodes3 = whiteNodes_roi_final - 1.5 * (grayNodes_roi_final - whiteNodes_roi_final);

        % regionPoints = [grayNodes_roi_final; interpolatedNodes1; interpolatedNodes2];
        regionPoints = [grayNodes_roi_final; whiteNodes_roi_final;interpolatedNodes1; interpolatedNodes2; interpolatedNodes3];

        % alpha shape -> boundary facets -> inpolyhedron
        shp = alphaShape(regionPoints, alphaRadius);
        [Connectivity_roi, vertices_roi] = boundaryFacets(shp);

        if isempty(Connectivity_roi) || isempty(vertices_roi)
            continue;
        end

        insideMask = inpolyhedron(Connectivity_roi, vertices_roi, brainPoints);
        insideMask_sum = insideMask_sum | insideMask;
    end

    % -------- 4) Apply final theta_y rules --------
    theta_y(unique(whiteElements(:))) = 0.0909;
    theta_y(insideMask_sum) = -0.4 * theta_y_temp(insideMask_sum);

    % output
    out.theta_y = theta_y;
    out.theta_y_temp = theta_y_temp;
    out.insideMask_sum = insideMask_sum;
    out.distToGraySurface = distToGraySurface;
end