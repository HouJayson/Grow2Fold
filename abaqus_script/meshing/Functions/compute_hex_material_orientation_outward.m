function [materialOrientations, materialOrientations_updated] = compute_hex_material_orientation_outward( ...
        pvertice_gray, pface_gray, pvertice_white, ...
        grayElements, brainElements, ...
        elemCentroids, elementRegions)

%--------------------------------------------------------------------------
% Compute outward material orientations for HEX8 brain elements
%
% Output:
%   materialOrientations_updated : (Nb x 9)
%       Row-wise [a1 a2 n] material frames for brainElements
%
% Notes:
%   - Normals are derived from the pial (gray) triangle surface
%   - Outward direction is enforced using white → gray direction
%   - elementRegions is defined per grayElements
%   - Regions 0 and 19 are assigned identity orientation
%--------------------------------------------------------------------------

Vg = pvertice_gray;
Fg = pface_gray;
Vw = pvertice_white;

% Ensure 1-based indexing for faces
if min(Fg(:)) == 0
    Fg = Fg + 1;
end

Ng = size(grayElements,1);
NvG = size(Vg,1);

% 1) Area-weighted vertex normals on gray surface
v1 = Vg(Fg(:,1),:);
v2 = Vg(Fg(:,2),:);
v3 = Vg(Fg(:,3),:);

fn = cross(v2 - v1, v3 - v1, 2);

I = [Fg(:,1); Fg(:,2); Fg(:,3)];
vn = zeros(NvG,3,'like',Vg);
for d = 1:3
    vn(:,d) = accumarray(I, repmat(fn(:,d),3,1), [NvG 1], @sum, 0);
end

vn_norm = vecnorm(vn,2,2);
vn(vn_norm>0,:) = vn(vn_norm>0,:) ./ vn_norm(vn_norm>0);

% 2) Enforce outward orientation (white → gray)
idxW = knnsearch(Vw, Vg, 'K', 1);
rad = Vg - Vw(idxW,:);
rad = rad ./ vecnorm(rad,2,2);

flipV = sum(vn .* rad, 2) < 0;
vn(flipV,:) = -vn(flipV,:);

% 3) Assign element normals (nearest gray vertex)
idxG = knnsearch(Vg, elemCentroids, 'K', 1);
nElem = vn(idxG,:);
nElem = nElem ./ vecnorm(nElem,2,2);

% 4) Build material frames (for ALL elements; no region exceptions)
materialOrientations = zeros(Ng,9,'like',Vg);

idx = (1:Ng).';              % <-- all gray elements
n = nElem(idx,:);

a1 = repmat([1 0 0], Ng, 1);
bad = abs(sum(a1.*n,2)) > 0.9;
a1(bad,:) = repmat([0 1 0], nnz(bad), 1);

a1 = a1 - (sum(a1.*n,2).*n);
a1 = a1 ./ vecnorm(a1,2,2);

a2 = cross(n, a1, 2);
a2 = a2 ./ vecnorm(a2,2,2);

materialOrientations(:,:) = [a1 a2 n];

%  4) Build material frames except for region 0 and 19
% materialOrientations = zeros(Ng,9,'like',Vg);
% 
% isGlobal = (elementRegions == 0) | (elementRegions == 19);
% materialOrientations(isGlobal,:) = repmat([1 0 0 0 1 0 0 0 1], nnz(isGlobal),1);
% 
% idx = find(~isGlobal);
% if ~isempty(idx)
%     n = nElem(idx,:);
%     a1 = repmat([1 0 0], numel(idx),1);
%     bad = abs(sum(a1.*n,2)) > 0.9;
%     a1(bad,:) = repmat([0 1 0], nnz(bad),1);
% 
%     a1 = a1 - (sum(a1.*n,2).*n);
%     a1 = a1 ./ vecnorm(a1,2,2);
% 
%     a2 = cross(n,a1,2);
%     a2 = a2 ./ vecnorm(a2,2,2);
% 
%     materialOrientations(idx,:) = [a1 a2 n];
% end

% 5) Map into brainElements
materialOrientations_updated = zeros(size(brainElements,1),9,'like',Vg);
materialOrientations_updated(:,[1 5 9]) = 1;

[~, loc] = ismember(grayElements, brainElements, 'rows');
valid = loc > 0;

materialOrientations_updated(loc(valid),:) = materialOrientations(valid,:);
end
