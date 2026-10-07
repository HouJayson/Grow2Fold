function [selectedNodeIndices, selectedElementIndices] = ...
    createWhiteMatterBoundarySet_face(brainPoints, whiteMatterElements, surfV, surfF, partialSelect)
% Select WM boundary nodes/elements by proximity to an interior part of a surface.
%
% Rules:
%   - Surface interior vertices: geodesic distance-to-boundary > boundaryMargin (along mesh)
%   - WM boundary nodes: distance to interior surface vertices <= distThreshold
%   - WM boundary elements: (CURRENT CODE) any node is boundary node
%       If you want ALL nodes boundary, switch "any" -> "all" below.
%
% partialSelect:
%   false (default): return all selected elements
%   true           : pick 4 elements near the centroid of selected elements

    if nargin < 5 || isempty(partialSelect)
        partialSelect = false;
    end

    distThreshold  = 0.2;
    boundaryMargin = 5.0;

    % ---- 1) Surface boundary vertices (open boundary edges) ----
    F = surfF;
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E, 2);
    [Eu, ~, ic] = unique(E, 'rows');
    edgeCount = accumarray(ic, 1);

    boundaryEdges = Eu(edgeCount == 1, :);
    boundaryVerts = unique(boundaryEdges(:));

    % If surface is closed (no boundary), treat all vertices as "interior"
    if isempty(boundaryVerts)
        interiorSurfV = surfV;
    else
        % ---- 2) Geodesic distance to boundary (along surface mesh) ----
        dGeo = geodesicDistanceToSources(surfV, surfF, boundaryVerts);

        interiorMask = dGeo > boundaryMargin;
        interiorSurfV = surfV(interiorMask, :);

        % Safety fallback
        if isempty(interiorSurfV)
            interiorSurfV = surfV;
        end
    end

    % Optional debug
    % scatter3(interiorSurfV(:,1),interiorSurfV(:,2),interiorSurfV(:,3), 5, 'filled');

    % ---- 3) WM nodes near interior surface vertices ----
    uniquePointIndices = unique(whiteMatterElements(:));
    wmPts = brainPoints(uniquePointIndices, :);

    [~, d] = knnsearch(interiorSurfV, wmPts);
    selectedNodeIndices = uniquePointIndices(d <= distThreshold);

    % ---- 4) Elements: select elements based on boundary-node membership ----
    isBoundaryNodeGlobal = false(size(brainPoints,1), 1);
    isBoundaryNodeGlobal(selectedNodeIndices) = true;

    % CURRENT behavior (your code): any node boundary => element selected
    selectedElementIndices = find(any(isBoundaryNodeGlobal(whiteMatterElements), 2));

    % If you truly want ALL nodes boundary => element selected, use:
    % selectedElementIndices = find(all(isBoundaryNodeGlobal(whiteMatterElements), 2));

    % ---- 5) Optional partial selection: pick 4 elements near center ----
    if partialSelect
        if numel(selectedElementIndices) <= 4
            % keep as is
        else
            elems = whiteMatterElements(selectedElementIndices, :);      % (nSel x nen)
            elemCentroids = squeeze(mean(reshape(brainPoints(elems(:),:), size(elems,1), size(elems,2), 3), 2));

            center = mean(elemCentroids, 1);
            dc = sqrt(sum((elemCentroids - center).^2, 2));

            [~, ord] = sort(dc, 'ascend');
            selectedElementIndices = selectedElementIndices(ord(1:1));
        end

        % recompute selectedNodeIndices to be consistent with 4 elements
        selectedNodeIndices = unique(whiteMatterElements(selectedElementIndices, :));
    end
end


function dMin = geodesicDistanceToSources(V, F, sourceVerts)
% Geodesic distance on a triangle mesh from multiple source vertices.
% Uses a weighted edge graph (weights = edge lengths) + shortest paths.

    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E, 2);
    E = unique(E, 'rows');

    w = sqrt(sum((V(E(:,1),:) - V(E(:,2),:)).^2, 2));
    G = graph(E(:,1), E(:,2), w, size(V,1));

    D = distances(G, sourceVerts);
    dMin = min(D, [], 1).';
end