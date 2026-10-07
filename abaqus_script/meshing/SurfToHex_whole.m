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

clear; clc;
currentDir = pwd;
addpath('Functions')

case_idx = 'subject';
sessionId_min='*';
ga_min=21;

fPath_storage = fullfile(currentDir, 'dataStorage',case_idx);
fPath_temp = fullfile(currentDir, 'temp',case_idx);

fPath_input = fullfile(currentDir, 'inputFiles','longitudinal','GW21',case_idx);

if ~exist(fPath_storage, 'dir') || ~exist(fPath_temp, 'dir')
    mkdir(fPath_storage);
    mkdir(fPath_temp);
end

fprintf('Processing sub-%s_session-%s\n',case_idx,sessionId_min);

casename=sprintf('%s.GW%d',case_idx,ga_min);
connectivityName = sprintf('Sconn_%s',case_idx);

Gray_surf = mvtk_read(fullfile(fPath_input,sprintf('sub-%s_ses-%s.GA%d.W.gray_hull_label.vtk',case_idx,sessionId_min,ga_min)));
White_surf = mvtk_read(fullfile(fPath_input,sprintf('sub-%s_ses-%s.GA%d.W.white_hull_reg_label.vtk',case_idx,sessionId_min,ga_min)));
inter_surf = mvtk_read(fullfile(fPath_input,sprintf('sub-%s_ses-%s.GA%d.slice.vtk',case_idx,sessionId_min,ga_min)));
skull_surf = mvtk_read(fullfile(fPath_input,sprintf('sub-%s_ses-%s.GA%d.inner_skull.vtk',case_idx,sessionId_min,ga_min)));

pvertice_gray = Gray_surf.vertices;
pvertice_white = White_surf.vertices;
pvertice_rigid = inter_surf.vertices;
pvertice_skull = skull_surf.vertices;

pface_gray = Gray_surf.faces;
pface_white = White_surf.faces;
pface_rigid = inter_surf.faces;
pface_skull = skull_surf.faces;

par_Clustering_gray = Gray_surf.par_Clustering;
par_Clustering_white = White_surf.par_Clustering;
par_huang_gray = Gray_surf.par_huang;
par_huang_white = White_surf.par_huang;
par_FS2009 = White_surf.par_FS2009;
par_MMP = White_surf.par_MMP;
%
%---------------------------------------------------------------------
%--------------------- Read mesh info from inp files -----------------
%---------------------------------------------------------------------
disp('Read mesh info from inp files')
fname = fullfile(fPath_input,sprintf('Brain_%s_Hex.inp',case_idx));
fid = fopen(fname, 'r');
raw_text = fread(fid, '*char')';
fclose(fid);
node_block = regexp(raw_text, '\*NODE[^\n]*\n(.*?)(?=\n\*)', 'tokens', 'once');
node_lines = strtrim(node_block{1});
node_data = textscan(node_lines, '%f%f%f%f', 'Delimiter', ',', 'CollectOutput', true);
nodes = node_data{1}; 
element_blocks = regexp(raw_text, '\*ELEMENT[^,\n]*, *TYPE=\w+, *ELSET=(\w+)[^\n]*\n(.*?)(?=(\n\*|$))','tokens');
elements = struct();
elset_names = cell(numel(element_blocks), 1);
for i = 1:numel(element_blocks)
    elset_name = element_blocks{i}{1};
    elset_names{i} = elset_name;
    raw_block = strtrim(element_blocks{i}{2});
    lines = splitlines(raw_block);
    joined = strjoin(lines, '\n');
    data = textscan(joined, '%d%d%d%d%d%d%d%d%d', 'Delimiter', ',', 'CollectOutput', true);
    elements.(elset_name) = data{1};  
end

old_ids = nodes(:,1); 
new_ids = (1:length(old_ids))';
nodes(:,1) = new_ids;
id_map = containers.Map(old_ids, new_ids);
for i = 1:numel(elset_names)
    name = elset_names{i};
    elems = elements.(name);
    elems(:,2:end) = arrayfun(@(id) id_map(id), elems(:,2:end));
    elements.(name) = elems;
end
brainPoints = nodes(:,2:4);
if isfield(elements, 'EB2')
    grayElements = elements.EB2(:, 2:9);
elseif isfield(elements, 'EB3')
    grayElements = elements.EB3(:, 2:9);
elseif isfield(elements, 'EB4')
    grayElements = elements.EB4(:, 2:9);
else
    error('Neither EB2 nor EB3 fields exist in the elements structure.');
end 
whiteElements = elements.EB1(:,2:9);
brainElements = [grayElements; whiteElements];
%%
%---------------------------------------------------------------------
%----------------------- Assign regional labels ----------------------
%---------------------------------------------------------------------
disp('Assign regional labels')

uniqueRegions = unique(par_Clustering_gray);
numRegions = length(uniqueRegions);
regionPoints = cell(numRegions, 1);
regionFaces  = cell(numRegions, 1);
minBound     = zeros(numRegions, 3);
maxBound     = zeros(numRegions, 3);

alphaVal = 1.5;   % adjust this value; smaller = tighter shape
for regionIdx = 1:numRegions
    regionLabel = uniqueRegions(regionIdx);
    grayNodesInRegion  = pvertice_gray(par_Clustering_gray   == regionLabel, :);
    whiteNodesInRegion = pvertice_white(par_Clustering_white == regionLabel, :);
    interpolatedNodes1 = 0.5 * (grayNodesInRegion + whiteNodesInRegion);
    pts = [grayNodesInRegion; whiteNodesInRegion;interpolatedNodes1];
    regionPoints{regionIdx} = pts;
    shp = alphaShape(pts(:,1), pts(:,2), pts(:,3), alphaVal);
    [faces, verts] = boundaryFacets(shp);
    regionFaces{regionIdx}  = faces;
    regionPoints{regionIdx} = verts;
    minBound(regionIdx,:) = min(verts, [], 1);
    maxBound(regionIdx,:) = max(verts, [], 1);
end

numElements = size(grayElements, 1);
nNodes = size(grayElements, 2);
elemNodeCoords = reshape(brainPoints(grayElements(:), :),numElements, nNodes, 3);
elemCentroids = squeeze(mean(elemNodeCoords, 2));
elementRegions = inf(numElements, 1);

for regionIdx = 1:numRegions
    faces = regionFaces{regionIdx};
    verts = regionPoints{regionIdx};
    rLab  = double(uniqueRegions(regionIdx));
    cand = all(elemCentroids >= minBound(regionIdx,:) & elemCentroids <= maxBound(regionIdx,:), 2);
    if ~any(cand)
        continue;
    end
    candIdx = find(cand);
    coords = reshape(elemNodeCoords(candIdx, :, :), [], 3);
    in = inpolyhedron(faces, verts, coords);
    in = reshape(in, numel(candIdx), nNodes);
    hitIdx = candIdx(any(in, 2));
    elementRegions(hitIdx) = min(elementRegions(hitIdx), rLab);
end
elementRegions(~isfinite(elementRegions)) = 0;


%---------------------------------------------------------------------
%------------- QC the assigned regional labels -----------------------
%---------------------------------------------------------------------
disp('QC assigned labels by neighbor majority vote (keep largest component)')
elementRegions = qc_labels_majority_vote(grayElements, elementRegions);
elementRegions_updated = zeros(length(brainElements),1) -1;
elementRegions_updated(ismember(brainElements,grayElements, 'rows'),:)=elementRegions;
disp('Region assignment completed.');

% filename=fullfile(fPath_temp, sprintf('%s_gray_regions.vtk',case_idx));
% writeVTKHexahedralMesh(filename, brainPoints, grayElements,elementRegions);


%---------------------------------------------------------------------
%---------- Calculate material orientation for each element ----------
%---------------------------------------------------------------------
disp('Calculate material orientation')
[materialOrientations, materialOrientations_updated] = compute_hex_material_orientation_outward(pvertice_gray, pface_gray, pvertice_white,grayElements, brainElements,elemCentroids, elementRegions);
disp('material orientations assignment completed.');


%----------------------------------------------------------------------------------------------
%---------- Create mesh and node set for ROIs (gray matter and white matter surface) ----------
%----------------------------------------------------------------------------------------------
disp('create mesh or node sets for boundary condition')
[boundaryNodeset, boundaryElementset] = createWhiteMatterBoundarySet_face(brainPoints, whiteElements,pvertice_rigid,pface_rigid,true);
nodeSets = cell(length(uniqueRegions), 1);
elementSets = cell(length(uniqueRegions), 1);
isPresentInGray = ismember(brainElements,grayElements, 'rows');
matchingRowIndices_gray = find(isPresentInGray);

for regionIdx = 0:numRegions-1
    elementIndices = find(elementRegions == regionIdx);
    regionElementNodes = grayElements(elementIndices, :);
    regionNodeIndices = unique(regionElementNodes(:));
    elementSets{regionIdx + 1} = matchingRowIndices_gray(elementIndices);
    nodeSets{regionIdx + 1} = regionNodeIndices;
end

elementPoints = brainPoints(brainElements, :);
elementSets{numRegions + 1}=sort(matchingRowIndices_gray);
nodeSets{numRegions + 1} = sort(unique(grayElements(:)));
isPresentInWhite = ismember(brainElements,whiteElements, 'rows');
matchingRowIndices_white = find(isPresentInWhite);
elementSets{numRegions + 2}=sort(matchingRowIndices_white);
nodeSets{numRegions + 2} = sort(unique(whiteElements(:)));
elementSets{numRegions + 3} =  matchingRowIndices_white(boundaryElementset);
nodeSets{numRegions + 3} = boundaryNodeset;

disp('create mesh or node sets for outermost surface of gray and white matter')
[surfaceFaces_gray, surfaceNodes_gray, surfaceElements_gray] = extractHexBoundaryFaces_fast(brainElements, brainPoints, pvertice_gray);
[surfaceFaces_white, surfaceNodes_white, surfaceElements_white] = extractHexBoundaryFaces_fast(whiteElements, brainPoints, pvertice_white);
[surfaceElement_gray_l, surfaceElement_gray_r,surfaceNodes_gray_l,surfaceNodes_gray_r, p_rigid] = splitSurfaceElementsByVTKPlane(brainPoints, brainElements, surfaceElements_gray, pvertice_rigid, pface_rigid, 1);
[rigid_nodes, rigid_elem, meta] = buildRigidWallBox(p_rigid);

nodeSets{numRegions + 4}=surfaceNodes_gray;
elementSets{numRegions + 4}=surfaceElements_gray;
nodeSets{numRegions + 5}=surfaceNodes_white;
elementSets{numRegions + 5}=matchingRowIndices_white(surfaceElements_white);
nodeSets{numRegions + 6}=rigid_nodes;
elementSets{numRegions + 6}=rigid_elem;
nodeSets{numRegions + 7}=[(1:size(pvertice_skull, 1))', pvertice_skull];
elementSets{numRegions + 7} = [(1:size(pface_skull, 1))', pface_skull];
nodeSets{numRegions + 8}=[p_rigid(1,:);mean(pvertice_skull,1)];

%%
%-----------------------------------------------------------------------------
%---------- Define the indicator function to apply the temperature -----------
%-----------------------------------------------------------------------------
%
disp('Define the indicator function to apply the temperature')
iterations = [0,6,3,2,2,5,4,   0,5,2,2,2,5,4]; % CC01019XX13

out = build_theta_indicator_temperature(brainPoints, whiteElements, pvertice_gray, pface_gray, pvertice_white, pvertice_gray,par_MMP, par_FS2009, iterations);
theta_y = out.theta_y;
nodeSets{numRegions + 9}=theta_y;
%%
disp('write files');
filename=fullfile(fPath_temp, sprintf('%s_temp_check.vtk',case_idx));
writeVTKHexahedralMesh(filename, brainPoints, brainElements,elementRegions_updated,materialOrientations_updated, theta_y);

disp('export the outer surface connectivity')
connectivityfname = fullfile(fPath_storage,[connectivityName,'.xlsx']);
exportSurfaceConnectivity_smooth_whole(surfaceNodes_gray, surfaceFaces_gray,surfaceNodes_white, surfaceFaces_white, connectivityfname);

disp('write abaqus inp file')
inpfname = fullfile(fPath_storage,[casename,'.par_Clustering.Hex_thk14_skull_2steps_stress.inp']);  
writeAbaqusInp_Hex_whole_rigbox_skull_2steps(inpfname, brainPoints, brainElements, materialOrientations, nodeSets, elementSets, numRegions)
