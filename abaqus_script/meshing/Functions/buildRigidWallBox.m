function [nodes, quads, meta] = buildRigidWallBox(p_rigid, planeSize, nGrid, halfThick)
% buildRigidWallBox
% Construct a "hollow box" rigid wall made of two parallel rectangular surfaces
% (front/back), each meshed into nGrid x nGrid quad elements.
%
% INPUTS
%   p_rigid   : (2x3)  p_rigid(1,:) = center point on plane
%                      p_rigid(2,:) = plane normal (does NOT need to be unit)
%   planeSize : (1x2)  [Lx Ly] size of rectangle in plane coordinates (default [200 200])
%   nGrid     : (1x2)  [nx ny] number of subdivisions (default [20 20])
%                      => nodes are (nx+1)x(ny+1) per side
%   halfThick : scalar half thickness along normal (default 0.1)
%               => surfaces at +/- halfThick along normal
%
% OUTPUTS
%   nodes : (N x 4)  [nodeID, x, y, z]
%   quads : (E x 5)  [elemID, n1, n2, n3, n4]  (4-node shell/quad connectivity)
%                    Ordering is consistent within each side; the back side is
%                    reversed to keep outward normals opposite.
%   meta  : struct with useful indexing maps
%
% EXAMPLE
%   p_rigid = [0 0 100; 1 0 0];  % center, normal
%   [nodes, quads] = buildRigidWallBox(p_rigid);
%   fprintf('%d nodes, %d quads\n', size(nodes,1), size(quads,1));

    if nargin < 2 || isempty(planeSize), planeSize = [200 200]; end
    if nargin < 3 || isempty(nGrid),     nGrid     = [20 20];   end
    if nargin < 4 || isempty(halfThick), halfThick = 0.15;       end

    c  = double(p_rigid(1,:));
    n0 = double(p_rigid(2,:));
    nNorm = norm(n0);
    if nNorm < 1e-12
        error('p_rigid(2,:) normal vector has near-zero magnitude.');
    end
    n = n0 / nNorm;

    % Build an orthonormal in-plane basis (u,v) for the plane
    % Pick a reference vector not parallel to n
    ref = [1 0 0];
    if abs(dot(ref,n)) > 0.9
        ref = [0 1 0];
    end
    u = cross(n, ref); u = u / norm(u);
    v = cross(n, u);   v = v / norm(v);

    % Grid parameters
    Lx = planeSize(1); Ly = planeSize(2);
    nx = nGrid(1);     ny = nGrid(2);

    xLin = linspace(-Lx/2, Lx/2, nx+1);   % along u
    yLin = linspace(-Ly/2, Ly/2, ny+1);   % along v
    [XX, YY] = meshgrid(xLin, yLin);      % (ny+1) x (nx+1)

    % Two surfaces: back (-) and front (+) along normal
    P0 = c + XX(:)*u + YY(:)*v;           % base plane points (N2D x 3)
    Pm = P0 - halfThick*n;                % back surface
    Pp = P0 + halfThick*n;                % front surface

    N2D = numel(XX);
    % Node IDs: first all back nodes, then all front nodes
    idBack  = (1:N2D).';
    idFront = (N2D+1:2*N2D).';

    nodesXYZ = [Pm; Pp];                  % (2*N2D x 3)
    nodes = [(1:2*N2D).', nodesXYZ];

    % Helper to convert (iy,ix) -> linear index in the meshgrid ordering
    % meshgrid with (y,x): linear index is sub2ind([ny+1, nx+1], iy, ix)
    idx2lin = @(iy, ix) sub2ind([ny+1, nx+1], iy, ix);

    % Build quad connectivity on one surface (back or front)
    % We’ll create quads with nodes (lower-left, lower-right, upper-right, upper-left)
    % in the (u,v) parametric space.
    E2D = nx * ny;
    quads_back  = zeros(E2D, 4);
    quads_front = zeros(E2D, 4);

    e = 0;
    for j = 1:ny
        for i = 1:nx
            e = e + 1;

            % corners in (j,i) cell
            ll = idx2lin(j,   i);
            lr = idx2lin(j,   i+1);
            ur = idx2lin(j+1, i+1);
            ul = idx2lin(j+1, i);

            % Back side uses same ordering
            quads_back(e,:) = idBack([ll lr ur ul]);

            % Front side reversed ordering so its normal points opposite (i.e., outward)
            % If you want both sides with same local normal, remove this reversal.
            quads_front(e,:) = idFront([ll ul ur lr]);
        end
    end

    % Combine element tables: first back, then front
    quads = [(1:(2*E2D)).', [quads_back; quads_front]];

    % Meta for convenience
    meta = struct();
    meta.n = n; meta.u = u; meta.v = v;
    meta.N2D = N2D; meta.E2D = E2D;
    meta.idBack  = reshape(idBack,  ny+1, nx+1);
    meta.idFront = reshape(idFront, ny+1, nx+1);
end
