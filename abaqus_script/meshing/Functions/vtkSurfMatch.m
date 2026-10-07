function [V_gray_reordered,F_updated] = vtkSurfMatch(V_gray,F_gray,V_refer, F_refer)

% ------ Surface matching (two vtk surface with the same face connectivity)-----
% it can also solve the issue that with a given surface connectivity, how to obtain the vertex ordering while keeping the surface morphology unchanged?
% -------(1 surface, two sets of face connectivity and vertex indexing)-----

disp('Obtaining the rough match based on the closet distance')
% Step 1: obtain the rough match based on the closet distance with the given referred surface 
N = size(V_refer, 1);
K = 5;  

gray_tree = KDTreeSearcher(V_gray);
[idx_all, ~] = knnsearch(gray_tree, V_refer, 'K', K); 
used_gray = false(N, 1);
idx_map = zeros(N, 1);

for i = 1:N
    for j = 1:K
        idx_candidate = idx_all(i, j);
        if ~used_gray(idx_candidate)
            idx_map(i) = idx_candidate;
            used_gray(idx_candidate) = true;
            break;
        end
    end
end

unmatched_refer = find(idx_map == 0);
unused_gray = find(~used_gray);

for i = 1:length(unmatched_refer)
    w_idx = unmatched_refer(i);
    
    distances = vecnorm(V_gray(unused_gray, :) - V_refer(w_idx, :), 2, 2);
    [~, min_idx] = min(distances);
    
    matched_idx = unused_gray(min_idx);
    idx_map(w_idx) = matched_idx;
    used_gray(matched_idx) = true;
    
    unused_gray(min_idx) = [];
end

V_gray_reordered = V_gray(idx_map, :);

disp('modifying incorrect connectivity in regional surface')
% Step 2: modifying incorrect connectivity in regional surface
tol = 1e-8;
N = size(V_gray, 1);
neighbors_ref = cell(N, 1);
for i = 1:size(F_refer, 1)
    face = F_refer(i, :);
    for j = 1:3
        v = face(j);
        others = face([1:j-1, j+1:end]);
        neighbors_ref{v} = [neighbors_ref{v}, others];
    end
end
for i = 1:N
    neighbors_ref{i} = unique(neighbors_ref{i});
end

neighbors_orig = cell(N, 1);
for i = 1:size(F_gray, 1)
    face = F_gray(i, :);
    for j = 1:3
        v = face(j);
        others = face([1:j-1, j+1:end]);
        neighbors_orig{v} = [neighbors_orig{v}, others];
    end
end
for i = 1:N
    neighbors_orig{i} = unique(neighbors_orig{i});
end

change_flag = true;
enable_multi_swap = true;

while change_flag
    change_flag = false;
    count = 0;
    V_gray_reordered0 = V_gray_reordered;

    for i = 1:N
        mapped_idx = find(all(abs(V_gray - V_gray_reordered(i,:)) < tol, 2), 1);
        if mapped_idx == 0
            continue;
        end

        nbr_ref_ids = neighbors_ref{i};
        nbr_mapped_ids = neighbors_orig{mapped_idx};

        [ref_coords_sorted, ref_sort_idx] = sortrows(V_gray(nbr_mapped_ids, :));
        [mapped_coords_sorted, mapped_sort_idx] = sortrows(V_gray_reordered0(nbr_ref_ids, :));

        if length(ref_sort_idx) ~= length(mapped_sort_idx)
            continue;
        end

        mismatch = any(abs(ref_coords_sorted - mapped_coords_sorted) > tol, 2);
        num_mismatch = sum(mismatch);

        % Adaptive behavior based on initial iteration count
        if enable_multi_swap && (num_mismatch == 2 || num_mismatch == 3)
            for j = find(mismatch)'
                vj = nbr_ref_ids(mapped_sort_idx(j));
                corrected_idx = nbr_mapped_ids(ref_sort_idx(j));
                match_coord = V_gray(corrected_idx, :);
                swap_target = find(all(abs(V_gray_reordered - match_coord) < 1e-6, 2), 1);
                temp = V_gray_reordered(vj, :);
                V_gray_reordered(vj, :) = V_gray_reordered(swap_target, :);
                V_gray_reordered(swap_target, :) = temp;
                change_flag = true;
                count= count +1;
            end
        elseif num_mismatch == 1
            j = find(mismatch);
            vj = nbr_ref_ids(mapped_sort_idx(j));
            corrected_idx = nbr_mapped_ids(ref_sort_idx(j));
            match_coord = V_gray(corrected_idx, :);
            swap_target = find(all(abs(V_gray_reordered - match_coord) < 1e-6, 2), 1);
            temp = V_gray_reordered(vj, :);
            V_gray_reordered(vj, :) = V_gray_reordered(swap_target, :);
            V_gray_reordered(swap_target, :) = temp;
            change_flag = true;
            count= count +1;
        end
    end

    fprintf('Total swaps: %d\n', count);
    enable_multi_swap = count > 2000;
end

F_updated = F_refer;

figure;
subplot(1,2,1)
trisurf(F_gray, V_gray(:,1), V_gray(:,2), V_gray(:,3));
title('Reference Surface'); axis equal;

subplot(1,2,2)
trisurf(F_refer, V_gray_reordered(:,1), V_gray_reordered(:,2), V_gray_reordered(:,3));
title('Reconstructed Surface'); axis equal;
end