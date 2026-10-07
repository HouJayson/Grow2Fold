function [surfaceFaces, surfaceNodes, surfaceElements] = extractHexBoundaryFaces_fast(hexElements, nodeCoordinates, pvertice_gray)

   distThresh = 1.5;  % keep your original threshold

    nE = size(hexElements,1);

    % --- Local face definition for a hex (same as your original order) ---
    fdef = [1 2 3 4;
            5 6 7 8;
            1 2 6 5;
            2 3 7 6;
            3 4 8 7;
            4 1 5 8];

    % --- Build all faces (6*nE x 4) ---
    allFaces = zeros(6*nE, 4, 'like', hexElements);
    for fid = 1:6
        rows = fid:6:(6*nE);
        allFaces(rows,:) = hexElements(:, fdef(fid,:));
    end

    % --- Find boundary faces: appear once ignoring orientation ---
    sortedFaces = sort(allFaces, 2);
    [~, ~, ic] = unique(sortedFaces, 'rows');      % group id per face
    counts = accumarray(ic, 1);
    isBoundary = counts(ic) == 1;

    bFaces = allFaces(isBoundary, :);              % candidate boundary faces

    % --- Compute face centroids (vectorized) ---
    c = ( nodeCoordinates(bFaces(:,1),:) + nodeCoordinates(bFaces(:,2),:) + ...
          nodeCoordinates(bFaces(:,3),:) + nodeCoordinates(bFaces(:,4),:) ) / 4;

    % --- Filter by distance-to-gray (pvertice_gray is large) ---
    keep = false(size(c,1),1);

    % Fast path: KD-tree nearest neighbor (Stats & ML Toolbox)
    useKD = false;
    if exist('createns','file') == 2 && exist('knnsearch','file') == 2
        useKD = true;
    end

    if useKD
        % Build KD-tree once, then query all centroids
        ns = createns(pvertice_gray, 'NSMethod', 'kdtree');
        % knnsearch returns nearest neighbor distance
        [~, d] = knnsearch(ns, c);
        keep = d <= distThresh;
    else
        % Fallback: chunked min squared distance (memory-safe, no toolbox)
        nB = size(c,1);
        minDist2 = inf(nB,1);
        c2 = sum(c.^2, 2);

        % Tune chunk size for your machine (larger = fewer loops but more RAM)
        chunk = 20000;

        for s = 1:chunk:size(pvertice_gray,1)
            e = min(s+chunk-1, size(pvertice_gray,1));
            G = pvertice_gray(s:e,:);

            g2 = sum(G.^2, 2)';            % 1 x m
            d2 = c2 + g2 - 2*(c*G');       % nB x m
            minDist2 = min(minDist2, min(d2, [], 2));
        end

        keep = minDist2 <= distThresh^2;
    end

    % --- Final surface faces ---
    surfaceFaces = bFaces(keep, :);

    % --- Surface nodes (fast logical mark) ---
    nNodes = size(nodeCoordinates,1);
    nodeMask = false(nNodes,1);
    nodeMask(surfaceFaces(:)) = true;
    surfaceNodes = find(nodeMask);

    % --- Surface elements (same logic as your original, but vectorized) ---
    % Original: element is surface if it contains >=4 surfaceNodes
    % Vectorized membership counts:
    isSurfNode = false(max(hexElements(:)), 1);
    isSurfNode(surfaceNodes) = true;

    surfCountPerElem = sum(isSurfNode(hexElements), 2);  % nE x 1
    surfaceElements = find(surfCountPerElem >= 4);
    surfaceElements = unique(surfaceElements);
end
