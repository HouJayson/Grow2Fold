function [surfaceElement_gray_l, surfaceElement_gray_r,surfaceNodes_gray_l, surfaceNodes_gray_r, P0123] = ...
    splitSurfaceElementsByVTKPlane( ...
    brainPoints, brainElements, surfaceElements_gray, vtkVertices, vtkFaces, dis_threshold)
    % splitSurfaceElementsByVTKPlane
    %
    % Outputs:
    %   P0123 : 5x3 CELL array (mixed numeric + formatted angle strings)
    %           Row 1: p0 = plane center (numeric)
    %           Row 2: p1 = p0 - n_k,  n_k = cross(nA, nB) (numeric)
    %           Row 3: p2 = [thetaDegStr, thetaDegStr, thetaDegStr] (string)
    %           Row 4: p3 = p0 - n_k_opposite, where n_k_opposite=cross(nA,-nB) (numeric)
    %           Row 5: p4 = [thetaOppDegStr, thetaOppDegStr, thetaOppDegStr] (string)
    %
    % Notes:
    %   nA is fixed as (0,1,0) (plane y=0, zx-plane)
    %   nB is the (unit) face normal of the VTK surface (best-fit oriented)
    
    %---------- sanity ----------
    if isempty(vtkFaces)
        error('vtkFaces must be provided to compute plane orientation robustly.');
    end
    
    P  = double(brainPoints);
    E  = double(brainElements);
    Vp = double(vtkVertices);
    Fp = double(vtkFaces);
    surfIdx = double(surfaceElements_gray(:));
    
    %---------- 1) centroids of gray surface elements ----------
    elemConn = E(surfIdx, :);
    elemConn(elemConn < 1) = NaN;
    
    nSurf = size(elemConn,1);
    centroids = zeros(nSurf,3);
    
    for i = 1:nSurf
        nodes = elemConn(i,:);
        nodes = nodes(~isnan(nodes));
        centroids(i,:) = mean(P(nodes,:),1);
    end
    
    %---------- 2) exclude elements close to vtk plane vertices ----------
    minDist = nearestPointDistance(centroids, Vp);
    keepMask = minDist >= dis_threshold;
    
    keptIdx = surfIdx(keepMask);
    keptCentroids = centroids(keepMask,:);
    
    %---------- 3) best-fit plane from vtk vertices ----------
    [planeCenter, planeNormal] = bestFitPlaneFromPoints(Vp, Fp); % oriented + unit
    p0 = planeCenter(:).';   % 1x3
    nB = planeNormal(:);     % 3x1
    nB = nB / norm(nB);
    

    %---------- 3.1) requested P0123 definition ---------- 
    nA = [0; 1; 0]; %plane A normal (y=0) 
    % axis direction (not normalized, per your instruction) 
    nk = cross(nA, nB); % n_k = nA x nB 
    p1 = (p0(:) - nk).'; % p1 = p0 - n_k 
    % rotation angle theta = atan2(||nA x nB||, nA · nB) 
    theta = atan2(norm(cross(nA, nB)), dot(nA, nB)); 
    thetaDeg = theta * 180/pi; 
    p2 = repmat(thetaDeg, 1, 3); 
    
    % opposite face normal 
    nBopp = -nB; 
    nkOpp = cross(nA, nBopp); % = cross(nA, -nB) = -nk 
    p3 = (p0(:) - nkOpp).'; % p3 = p0 - n_k_opposite (= p0 + nk) 
    thetaOpp = atan2(norm(cross(nA, nBopp)), dot(nA, nBopp)); 
    thetaOppDeg = thetaOpp * 180/pi; p4 = repmat(thetaOppDeg, 1, 3);

    P0123 = zeros(2,3);
    P0123(1,:) = p0;
    P0123(2,:) = nB';     
    % P0123(3,:) = p1;
    % P0123(4,:) = p2;
    % P0123(5,:) = p3;
    % P0123(6,:) = p4;

    
    %---------- 4) split left / right ----------
    if isempty(keptIdx)
        surfaceElement_gray_l  = [];
        surfaceElement_gray_r  = [];
    else
        signedDist = (keptCentroids - p0) * nB;  % signed distance to plane
        % keep your convention (left: >0, right: <0)
        surfaceElement_gray_l  = keptIdx(signedDist > 0);
        surfaceElement_gray_r  = keptIdx(signedDist < 0);
    end
    
    %---------- 5) surface nodes corresponding to left/right element sets ----------
    surfaceNodes_gray_l = elementsToUniqueNodes(E, surfaceElement_gray_l);
    surfaceNodes_gray_r = elementsToUniqueNodes(E, surfaceElement_gray_r);
    end
    
    % ================= helper functions =================
    
    function minDist = nearestPointDistance(Q, P)
    try
        [~, minDist] = knnsearch(P, Q);
    catch
        N = size(Q,1);
        minDist = zeros(N,1);
        chunk = 5000;
        for s = 1:chunk:N
            t = min(s+chunk-1, N);
            D = pdist2(Q(s:t,:), P);
            minDist(s:t) = min(D, [], 2);
        end
    end
    end
    
    function [center, normal] = bestFitPlaneFromPoints(pts, faces)
    center = mean(pts,1);
    X = pts - center;
    
    C = (X' * X) / size(X,1);
    [V,D] = eig(C);
    [~,idx] = min(diag(D));
    normal = V(:,idx);
    normal = normal / norm(normal);
    
    % orient consistently using face normals
    Fn = estimateMeshNormal(pts, faces);
    if dot(normal, Fn) < 0
        normal = -normal;
    end
    end
    
    function nAvg = estimateMeshNormal(V, F)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);
    
    fn = cross(v2 - v1, v3 - v1, 2);
    nSum = sum(fn,1);
    
    if norm(nSum) < 1e-12
        nAvg = [0;0;1];
    else
        nAvg = nSum(:) / norm(nSum);
    end
    end
    
    function nodeIdx = elementsToUniqueNodes(E, elemIdx)
    % Convert a list of element indices into unique global node indices.
    if isempty(elemIdx)
        nodeIdx = [];
        return;
    end
    conn = E(double(elemIdx), :);
    conn = conn(:);
    conn = conn(conn > 0);            % drop zeros if any
    conn = conn(~isnan(conn));        % drop NaNs if any
    nodeIdx = unique(conn);
    end
