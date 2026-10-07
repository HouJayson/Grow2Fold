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
% Refactored NMF/opNMF/spectral cortical parcellation pipeline
% Added:
%   - graph-regularized opNMF backend
%   - spectral clustering backend
%   - soft hierarchy prior for spectral clustering across k
%   - hierarchy score in model selection
%   - explicit sorting by increasing gestational age before forming V
clear; clc; close("all");
rng(1);

%% ==================== USER SETTINGS ====================
cfg = struct();

% ---------- Input ----------
cfg.data_root = "root-path";
cfg.pattern.L = "path\*.lh.midthickness.fsavg10k_smth.vtk";
cfg.pattern.R = "path\*.rh.midthickness.fsavg10k_smth.vtk";
cfg.inflated.L = "path\avg.GA40.lh.inflated.fsavg10k_new.vtk";
cfg.inflated.R = "path\avg.GA40.rh.inflated.fsavg10k_new.vtk";

% ---------- Subject / age filter ----------
cfg.min_GA = 21;
cfg.max_GA = 44;

% ---------- Preprocessing ----------
cfg.exclude_medialwall = true;
cfg.medial_name = "par_huang";  % par_huang
cfg.eps_floor = 1e-12;

cfg.normalize_columns = true;
cfg.column_norm_type = "mean";   % "none" | "l1" | "l2" | "mean" | "log_mean" | "sqrt_mean" | "quantile_mean" | "range01_col" | "zscore_pos" | "robust_zscore_pos"
cfg.quantile_q = 0.95;
cfg.zscore_shift = "min";        % "min" | "3sd"
cfg.robust_mad_scale = 1.4826;

cfg.smooth_enable = true;
cfg.smooth_iters = 5;

% ---------- Partition backend ----------
cfg.nmf_backend = "spectral";    % "opnmf" | "opnmf_graph" | "nnmf" | "spectral"

cfg.nRep_search = 1;
cfg.nRep_final  = 2;
cfg.maxIter_search = 50000;
cfg.maxIter_final  = 100000;
cfg.tolFun = 1e-5;
cfg.nnmf_algorithm = "mult";

% ---------- Spectral clustering ----------
cfg.spec_similarity = "traj_ref_fused";   % "traj_corr" | "traj_ref_fused"
cfg.spec_ref_n = 320;                     % reference vertices for ref-profile similarity
cfg.spec_knn = 0;                        % final KNN sparsification; [] or 0 = dense; additional graph sparsification after fused simularity matrix
cfg.spec_use_normalized_laplacian = true;
cfg.spec_kmeans_reps = 20;                % number of random initializations
cfg.spec_kmeans_maxIter = 2000;           % max iterations per k-means run

% correlation type for similarity computation
cfg.spec_corr_type = "pearson";           % "pearson" | "spearman" | "kendall"

% fusion type for traj_ref_fused
cfg.spec_fusion_method = "snf";       % "average" | "snf"
cfg.spec_fuse_alpha = 0.5;            % only used when spec_fusion_method = "average"
cfg.spec_snf_niter = 3;              % Number of SNF iterations
cfg.spec_snf_knn = 50;                % vertex keeps top ** neighbors during fusion

% check the similarity matrix
cfg.spec_save_similarity = true;
cfg.spec_save_similarity_dir = 'F:\NMF\temp';
cfg.spec_save_similarity_format = 'png';
cfg.spec_save_similarity_max_n = 20;
cfg.spec_save_reorder_by_label = true;
cfg.spec_save_snf_each_iter = false;
cfg.spec_save_snf_kernels = false;   % save P and Q

% ---------- Graph regularization ----------
cfg.graph_lambda = 0.1;                   % try: 0, 0.01, 0.05, 0.1, 0.5
cfg.graph_use_normalized_laplacian = false;
cfg.graph_check_obj_every = 100;
cfg.orth_beta = 0.05;   % try 0, 0.01, 0.05, 0.1, 0.5

% ---------- k selection ----------
cfg.k_range = 2:2;
cfg.k_min_select = 10;
cfg.use_shared_k = false;

% Metric 1: silhouette on cortex-only sampled vertices
cfg.sil_residualize = false;
cfg.sil_sample_n = 10242;  % 2562, 10242 , 32492
cfg.sil_distance = 'correlation';         % 'euclidean' | 'correlation'

% Metric 2: reconstruction error
cfg.recon_error_type = "col_rel_sse";     % "abs_sse" | "rel_fro" | "rel_sse" | "col_rel_sse"

% Metric 3: instability
cfg.stab_enable = true;
cfg.stab_mode = "paper";                  % "paper" | "split" | "none"
cfg.stab_metric = "ari";                  % "ari" | "nmi"
cfg.stab_nSplits = 5;
cfg.stab_min_scans = 100;
cfg.stab_weight = 0.5;

cfg.write_stability_vtk = false;
cfg.stability_vtk_mode = "all";           % "all" | "best" | "list"
cfg.stability_vtk_k_list = [];
cfg.stability_vtk_dirname = "vtk_stab";

% ---------- Label cleanup ----------
cfg.cleanup_enable = true;
cfg.island_cleanup_mode = "both";         % "component" | "majority" | "both"
cfg.min_island_size = 200;
cfg.majority_vote_iters = 10;
cfg.majority_vote_frac = 0.55;

% ---------- Output ----------
cfg.out_dir = fullfile(pwd, "nmf_out_snf");
cfg.write_vtk_labels = true;
cfg.vtk_mode = "all";                     % "best" | "all" | "list"
cfg.vtk_k_list = [];
cfg.vtk_sweep_dirname = "vtk_sweep";

addpath("F:\NMF\brainparts");
if ~exist(cfg.out_dir, 'dir'); mkdir(cfg.out_dir); end

%% ==================== RUN ====================
fprintf("\n========== LEFT (metrics) ==========\n");
resL = run_one_hemi("L", cfg, true, []);

% fprintf("\n========== RIGHT (metrics) ==========\n");
% resR = run_one_hemi("R", cfg, true, []);
% 
% if cfg.use_shared_k
%     shared_k = pick_shared_k(resL.metrics, resR.metrics, cfg);
% else
%     shared_k = [];
% end

% fprintf("\n========== FINAL FIT ==========\n");
% run_one_hemi("L", cfg, false, shared_k);
% run_one_hemi("R", cfg, false, shared_k);

% Define parameter ranges
% snf_knn_values = [50 500 1000];
% knn_values = 20;
% 
% % Outer loop for spec_ref_n
% for ref_n = snf_knn_values
%     cfg.spec_snf_knn = ref_n;
%     % Inner loop for spec_knn
%     for knn = knn_values
%         cfg.spec_knn = knn;
%         cfg.vtk_sweep_dirname = sprintf('vtk_sweep_200removal_snf_knn_values%d_knn%d',ref_n,knn);
% 
%         fprintf("\n========== LEFT (metrics): snf_knn_values=%d, spec_knn=%d ==========\n",ref_n,knn);
%         resL = run_one_hemi("L", cfg, true, []);
% 
%         fprintf("\n========== RIGHT (metrics): snf_knn_values=%d, spec_knn=%d ==========\n",ref_n,knn);
%         resR = run_one_hemi("R", cfg, true, []);
%     end
% end


% fprintf("\n========== FINAL FIT ==========\n");
% run_one_hemi("L", cfg, false, shared_k);
% run_one_hemi("R", cfg, false, shared_k);

% Define parameter ranges
% snf_niter_values = 20;
% knn_values = 25;
% 
% % Outer loop for spec_ref_n
% for ref_n = snf_niter_values
%     cfg.spec_snf_niter = ref_n;
%     % Inner loop for spec_knn
%     for knn = knn_values
%         cfg.spec_knn = knn;
%         cfg.vtk_sweep_dirname = sprintf('vtk_sweep_3k_snf_niter_values%d_knn%d',ref_n,knn);
% 
%         fprintf("\n========== LEFT (metrics): snf_niter_values=%d, spec_knn=%d ==========\n",ref_n,knn);
%         resL = run_one_hemi("L", cfg, true, []);
% 
%         fprintf("\n========== RIGHT (metrics): snf_niter_values=%d, spec_knn=%d ==========\n",ref_n,knn);
%         resR = run_one_hemi("R", cfg, true, []);
%     end
% end

fprintf("\nAll done.\n");

%% ==================== MAIN PIPELINE ====================
function out = run_one_hemi(hemi, cfg, metrics_only, forced_k)
    out = struct();
    out.hemi = hemi;

    out_dir = fullfile(cfg.out_dir, "hemi_" + hemi, cfg.vtk_sweep_dirname);
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    files = find_surface_files(cfg.data_root, cfg.pattern.(char(hemi)));
    files = filter_files_by_GA(files, cfg.min_GA, cfg.max_GA);
    if isempty(files)
        error("No scans remain for hemi %s after GA filtering.", hemi);
    end

    % ---- sort scans by increasing gestational age ----
    [files, GA_sorted] = sort_files_by_GA(files);
    fprintf('[Hemi %s] Files sorted by increasing GA. Range: %.2f -> %.2f (%d scans)\n', ...
        hemi, GA_sorted(1), GA_sorted(end), numel(GA_sorted));

    [V_full, mesh] = build_area_matrix(files, cfg);
    mask = mesh.mask_cortex;
    
    % Build graph on full mesh, then restrict to cortex-only analysis vertices
    A_full = build_sparse_adjacency(mesh.nV_full, mesh.faces);
    A_ctx = A_full(mask, mask);
    
    % cortex-only coordinates for spectral reference sampling
    vertex_coords = mesh.vertices(mask, :);

    % ---- DEBUG: visualize reference sampling ----
    if strcmpi(cfg.nmf_backend, "spectral")
        [V_ctx, F_ctx] = extract_cortex_mesh(mesh.vertices, mesh.faces, mask);
        
        ref_idx = uniform_ref_indices_fps(V_ctx, cfg.spec_ref_n, [], []);
        
        write_reference_points_vtk( ...
            fullfile(out_dir, sprintf("ref_points_hemi_%s.vtk", hemi)), ...
            V_ctx, F_ctx, ref_idx);
    end
    
    % Use cortex-only matrix for analysis; normalize after masking.
    V = preprocess_V(V_full(mask, :), cfg);

    fprintf("[Hemi %s] V_full = %d x %d, V_cortex = %d x %d\n", ...
        hemi, size(V_full,1), size(V_full,2), size(V,1), size(V,2));

    GA = arrayfun(@parse_GA_from_filename, files);
    fprintf('[Hemi %s] First 10 sorted GAs: ', hemi);
    fprintf('%.2f ', GA(1:min(10,end)));
    fprintf('\n');

    if cfg.sil_residualize
        V_sil = residualize_by_GA(V, GA);
    else
        V_sil = V;
    end

    if isempty(forced_k)
        [best_k, metrics] = select_k_joint(V, V_sil, A_ctx, vertex_coords, mesh, hemi, cfg, out_dir);
        plot_k_selection(metrics, hemi, out_dir);
    else
        best_k = forced_k;
        metrics = struct('best_k', best_k, 'k_range', cfg.k_range(:), 'k_min_select', cfg.k_min_select);
        plot_k_selection_forced(out_dir, hemi, best_k);
    end

    out.metrics = metrics;
    if metrics_only
        fprintf("[Hemi %s] Metrics-only mode.\n", hemi);
        return;
    end

    fprintf("[Hemi %s] Final %s fit with k = %d ...\n", hemi, upper(cfg.nmf_backend), best_k);

    [labels_ctx, model] = fit_partition_backend(V, A_ctx, vertex_coords, best_k, cfg, "final");

    labels_full = zeros(mesh.nV_full, 1, 'int32');
    labels_full(mask) = labels_ctx;

    if cfg.cleanup_enable
        labels_full = cleanup_labels_locked(labels_full, mesh.faces, mask, cfg);
    end

    save_outputs(model, labels_ctx, labels_full, best_k, metrics, cfg, files, mesh, hemi, out_dir);
end

%% ==================== K SELECTION ====================
function [best_k, metrics] = select_k_joint(V_nmf, V_sil, A_ctx, vertex_coords, mesh, hemi, cfg, out_dir)
    k_range = cfg.k_range(:);
    nK = numel(k_range);
    nV = size(V_nmf, 1);

    sample_n = min(cfg.sil_sample_n, nV);
    sample_idx = randperm(nV, sample_n);
    V_sil_samp = V_sil(sample_idx, :);

    recon = nan(nK, 1);
    sil = nan(nK, 1);
    instability = nan(nK, 1);

    write_sweep_vtk = should_write_sweep_vtk(cfg, mesh);
    if write_sweep_vtk
        vtk_dir = fullfile(out_dir, cfg.vtk_sweep_dirname);
        if ~exist(vtk_dir, 'dir'); mkdir(vtk_dir); end
    end

    for i = 1:nK
        k = k_range(i);
        fprintf("[Hemi %s] Testing k = %d (%d/%d) ...\n", hemi, k, i, nK);

        [labels_full_ctx, model] = fit_partition_backend(V_nmf, A_ctx, vertex_coords, k, cfg, "search");

        if isfield(model, 'W') && isfield(model, 'H') && ~isempty(model.W) && ~isempty(model.H)
            R = V_nmf - model.W * model.H;
            recon(i) = compute_recon_error(V_nmf, R, cfg.recon_error_type);
        else
            recon(i) = NaN;
        end

        labels_s = labels_full_ctx(sample_idx);
        sil(i) = compute_mean_silhouette(V_sil_samp, labels_s, cfg);

        if cfg.stab_enable && cfg.stab_mode ~= "none"
            switch lower(cfg.stab_mode)
                case "paper"
                    instability(i) = paper_style_instability(V_nmf, A_ctx, vertex_coords, k, cfg, labels_full_ctx, mesh, hemi, out_dir);
                case "split"
                    instability(i) = 1 - split_half_stability(V_nmf, A_ctx, vertex_coords, k, cfg);
                otherwise
                    error("Unknown cfg.stab_mode: %s", cfg.stab_mode);
            end
        end

        fprintf("    Recon = %.4e | Sil = %.4f | Instability = %.4f\n", ...
            recon(i), sil(i), instability(i));

        if write_sweep_vtk && should_write_this_k(k, cfg)
            labels_full = zeros(mesh.nV_full, 1, 'int32');
            labels_full(mesh.mask_cortex) = labels_full_ctx;
            if cfg.cleanup_enable
                labels_full = cleanup_labels_locked(labels_full, mesh.faces, mesh.mask_cortex, cfg);
            end
            [Vinf, Finf] = load_inflated_mesh(hemi, cfg, mesh);
            write_vtk_point_scalar(fullfile(vtk_dir, sprintf("labels_hemi_%s_k%02d_sweep.vtk", hemi, k)), ...
                Vinf, Finf, labels_full, "NMF_label");
        end
    end

    recon_z = zscore_nan_safe(recon);
    sil_z = zscore_nan_safe(sil);

    use_recon = ~strcmpi(cfg.nmf_backend, "spectral") && any(isfinite(recon));

    if cfg.stab_enable && cfg.stab_mode ~= "none"
        instab_z = zscore_nan_safe(instability);
        if use_recon
            joint = sil_z - recon_z - cfg.stab_weight * instab_z;
        else
            joint = sil_z - cfg.stab_weight * instab_z;
        end
    else
        if use_recon
            joint = sil_z - recon_z;
        else
            joint = sil_z;
        end
    end

    joint(k_range < cfg.k_min_select) = -Inf;
    [~, idx] = max(joint);
    best_k = k_range(idx);

    metrics = struct();
    metrics.k_range = k_range;
    metrics.reconstruction_error = recon;
    metrics.silhouette_mean = sil;
    metrics.instability = instability;
    metrics.recon_z = recon_z;
    metrics.sil_z = sil_z;
    metrics.instab_z = [];
    if cfg.stab_enable && cfg.stab_mode ~= "none"
        metrics.instab_z = instab_z;
    end
    metrics.joint_score = joint;
    metrics.sample_idx = sample_idx;
    metrics.best_k = best_k;
    metrics.k_min_select = cfg.k_min_select;

    save(fullfile(out_dir, "k_selection_metrics.mat"), "metrics");
    writetable(table(k_range, recon, sil, instability, joint, ...
        'VariableNames', {'k','reconstruction_error','silhouette_mean','instability','joint_score'}), ...
        fullfile(out_dir, "k_selection_summary.csv"));

    fprintf("[Hemi %s] Selected best k = %d (k >= %d)\n", hemi, best_k, cfg.k_min_select);
end

%% ==================== PARTITION BACKEND ====================
function [labels, model] = fit_partition_backend(V, A, vertex_coords, k, cfg, phase)
    switch lower(cfg.nmf_backend)
        case {"opnmf","opnmf_graph","nnmf"}
            [W, H] = fit_nmf_backend(V, A, k, cfg, phase);
            labels = hard_labels_from_W(W);
            model = struct('W', W, 'H', H);

        case "spectral"
            labels = spectral_cluster(V, vertex_coords, k, cfg);
            model = struct('W', [], 'H', []);

        otherwise
            error("Unknown cfg.nmf_backend: %s", cfg.nmf_backend);
    end
end

function [Wbest, Hbest] = fit_nmf_backend(V, A, k, cfg, phase)
    switch lower(string(phase))
        case "search"
            nRep = cfg.nRep_search;
            maxIter = cfg.maxIter_search;
        otherwise
            nRep = cfg.nRep_final;
            maxIter = cfg.maxIter_final;
    end

    switch lower(cfg.nmf_backend)
        case "opnmf"
            [Wbest, Hbest] = opnmf_replicates(V, k, nRep, maxIter, cfg.tolFun);

        case "opnmf_graph"
            [Wbest, Hbest] = opnmf_graph_replicates(V, A, k, nRep, maxIter, cfg.tolFun, cfg.graph_lambda, cfg.orth_beta, cfg);

        case "nnmf"
            [Wbest, Hbest] = nnmf_replicates(V, k, nRep, maxIter, cfg.tolFun, cfg.nnmf_algorithm);

        otherwise
            error("Unknown cfg.nmf_backend: %s", cfg.nmf_backend);
    end
end

function [Wbest, Hbest, bestErr] = opnmf_replicates(V, k, nRep, maxIter, tol)
    bestErr = inf;
    Wbest = [];
    Hbest = [];
    for r = 1:nRep
        [W, ~] = opnmf(V, k, [], 0, maxIter, tol);
        H = W' * V;
        err = norm(V - W * H, 'fro')^2;
        if err < bestErr
            bestErr = err;
            Wbest = W;
            Hbest = H;
        end
    end
end

function [Wbest, Hbest, bestErr] = opnmf_graph_replicates(V, A, k, nRep, maxIter, tol, lambda, beta_orth, cfg)
    bestErr = inf;
    Wbest = [];
    Hbest = [];
    for r = 1:nRep
        [W, ~, info] = opnmf_graph(V, k, A, lambda, beta_orth, [], 0, maxIter, tol, cfg);
        H = W' * V;

        if ~isempty(info.obj) && any(isfinite(info.obj))
            idx = find(isfinite(info.obj), 1, 'last');
            err = info.obj(idx);
        else
            err = norm(V - W * H, 'fro')^2;
        end

        if err < bestErr
            bestErr = err;
            Wbest = W;
            Hbest = H;
        end
    end
end

function [Wbest, Hbest] = nnmf_replicates(V, k, nRep, maxIter, tol, alg)
    V = max(V, 0);
    opts = statset('MaxIter', maxIter, 'TolFun', tol, 'Display', 'off');
    [Wbest, Hbest] = nnmf(V, k, 'algorithm', char(alg), 'replicates', nRep, 'options', opts);
end

%% ==================== DATA BUILDING ====================
function [files_out, ga_sorted] = sort_files_by_GA(files_in)
    ga = arrayfun(@parse_GA_from_filename, files_in);

    if any(isnan(ga))
        error('sort_files_by_GA: some files have unreadable GA after filtering.');
    end

    [ga_sorted, ord] = sort(ga, 'ascend');
    files_out = files_in(ord);
end

function files = find_surface_files(root_dir, pattern)
    d = dir(fullfile(char(root_dir), char(pattern)));
    files = strings(0,1);
    for i = 1:numel(d)
        if ~d(i).isdir
            files(end+1,1) = string(fullfile(d(i).folder, d(i).name)); %#ok<AGROW>
        end
    end
end

function [V_full, mesh] = build_area_matrix(scan_files, cfg)
    [V0, F0, S0] = read_surface_vtk(scan_files(1));
    nV = size(V0,1);

    mesh.vertices = V0;
    mesh.faces = F0;
    mesh.nV_full = nV;

    if cfg.exclude_medialwall
        if ~isfield(S0, cfg.medial_name)
            error("Field %s not found in VTK scalars.", cfg.medial_name);
        end
        mw = double(S0.(cfg.medial_name)(:));
        mask = (mw ~= 0);
        mesh.(cfg.medial_name) = mw;
    else
        mask = true(nV,1);
    end
    mesh.mask_cortex = mask;

    nScans = numel(scan_files);
    V_full = zeros(nV, nScans, 'double');
    V_full(:,1) = vertex_areas_from_mesh(V0, F0);

    for s = 2:nScans
        [Vs, Fs] = read_surface_vtk(scan_files(s));
        if size(Vs,1) ~= nV || ~isequal(Fs, F0)
            error("Topology mismatch at scan %d.", s);
        end
        V_full(:,s) = vertex_areas_from_mesh(Vs, Fs);
        if mod(s,100) == 0 || s == nScans
            fprintf("Loaded %d/%d scans\n", s, nScans);
        end
    end

    if cfg.smooth_enable
        V_full = smooth_vertex_maps_masked_full(V_full, F0, mask, cfg.smooth_iters);
    end
end

function V = preprocess_V(V, cfg)
    V = double(V);
    V(~isfinite(V)) = 0;
    V(V < 0) = 0;
    V(V < cfg.eps_floor) = 0;

    norm_type = "none";
    if isfield(cfg, 'normalize_columns') && cfg.normalize_columns
        norm_type = lower(string(cfg.column_norm_type));
    end

    switch norm_type
        case "none"
        case "l1"
            V = V ./ (sum(V,1) + cfg.eps_floor);
        case "l2"
            V = V ./ (sqrt(sum(V.^2,1)) + cfg.eps_floor);
        case "mean"
            V = V ./ (mean(V,1) + cfg.eps_floor);
        case "log_mean"
            V = log1p(V);
            V = V ./ (mean(V,1) + cfg.eps_floor);
        case "sqrt_mean"
            V = sqrt(V);
            V = V ./ (mean(V,1) + cfg.eps_floor);
        case "quantile_mean"
            q = 0.95;
            if isfield(cfg, 'quantile_q') && ~isempty(cfg.quantile_q)
                q = cfg.quantile_q;
            end
            denom = prctile(V, 100*q, 1);
            V = V ./ (denom + cfg.eps_floor);
            V = V ./ (mean(V,1) + cfg.eps_floor);
        case "range01_col"
            cmin = min(V, [], 1);
            cmax = max(V, [], 1);
            V = (V - cmin) ./ (cmax - cmin + cfg.eps_floor);
        case "zscore_pos"
            mu = mean(V, 1);
            sd = std(V, 0, 1);
            sd(sd < cfg.eps_floor) = 1;
            V = (V - mu) ./ sd;
            shift_mode = "min";
            if isfield(cfg, 'zscore_shift') && ~isempty(cfg.zscore_shift)
                shift_mode = lower(string(cfg.zscore_shift));
            end
            switch shift_mode
                case "min"
                    V = V - min(V, [], 1);
                case "3sd"
                    V = V + 3;
                    V(V < 0) = 0;
                otherwise
                    error("Unknown cfg.zscore_shift: %s", shift_mode);
            end
            V = V ./ (mean(V,1) + cfg.eps_floor);
        case "robust_zscore_pos"
            medv = median(V, 1);
            madv = median(abs(V - medv), 1);

            mad_scale = 1.4826;
            if isfield(cfg, 'robust_mad_scale') && ~isempty(cfg.robust_mad_scale)
                mad_scale = cfg.robust_mad_scale;
            end
            madv = mad_scale * madv;
            madv(madv < cfg.eps_floor) = 1;

            V = (V - medv) ./ madv;
            V = V - min(V, [], 1);
            V = V ./ (mean(V,1) + cfg.eps_floor);
        otherwise
            error("Unknown normalization type: %s", cfg.column_norm_type);
    end
    V(~isfinite(V)) = 0;
    V(V < 0) = 0;
end

function Vres = residualize_by_GA(V, GA)
    X = [ones(numel(GA),1), GA(:)];
    B = (X' * X) \ (X' * V');
    Vhat = (X * B)';
    Vres = V - Vhat;
end

%% ==================== MESH / VTK ====================
function [V, F, S] = read_surface_vtk(vtk_path)
    S = mvtk_read(char(vtk_path));

    if isfield(S, 'vertices'); V = double(S.vertices);
    elseif isfield(S, 'points'); V = double(S.points);
    else; error('No vertices/points field found.');
    end

    if isfield(S, 'faces'); F = double(S.faces);
    elseif isfield(S, 'triangles'); F = double(S.triangles);
    elseif isfield(S, 'polygons'); F = double(S.polygons);
    else; error('No faces/triangles/polygons field found.');
    end

    if min(F(:)) == 0; F = F + 1; end
    if size(F,2) ~= 3; error('Mesh must be triangular.'); end
end

function a = vertex_areas_from_mesh(V, F)
    v1 = V(F(:,1),:); v2 = V(F(:,2),:); v3 = V(F(:,3),:);
    triA = 0.5 * vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);
    nV = size(V,1);
    a = accumarray(F(:,1), triA/3, [nV,1], @sum, 0) + ...
        accumarray(F(:,2), triA/3, [nV,1], @sum, 0) + ...
        accumarray(F(:,3), triA/3, [nV,1], @sum, 0);
    a(~isfinite(a)) = 0;
    a(a < 0) = 0;
end

function [Vinf, Finf] = load_inflated_mesh(hemi, cfg, mesh)
    try
        [Vinf, Finf] = read_surface_vtk(fullfile(cfg.data_root, cfg.inflated.(char(hemi))));
        if ~isequal(Finf, mesh.faces)
            error('Inflated mesh topology mismatch.');
        end
    catch
        Vinf = mesh.vertices;
        Finf = mesh.faces;
    end
end

function A = build_sparse_adjacency(nV, F)
    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    A = sparse(E(:,1), E(:,2), 1, nV, nV);
    A = A + A';
    A(A > 0) = 1;
    A = A - diag(diag(A));
end

function tf = should_write_stability_this_k(k, cfg)
    switch lower(cfg.stability_vtk_mode)
        case 'all'
            tf = true;
        case 'best'
            tf = false;
        case 'list'
            tf = ismember(k, cfg.stability_vtk_k_list);
        otherwise
            tf = false;
    end
end

function write_vtk_point_scalar(out_path, V, F, scalar, scalar_name)
    fid = fopen(char(out_path), 'w');
    if fid < 0; error('Cannot open file for writing: %s', out_path); end

    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'NMF labels\nASCII\nDATASET POLYDATA\n');
    fprintf(fid, 'POINTS %d float\n', size(V,1));
    fprintf(fid, '%f %f %f\n', V');

    fprintf(fid, 'POLYGONS %d %d\n', size(F,1), size(F,1)*4);
    fprintf(fid, '3 %d %d %d\n', (F-1)');

    fprintf(fid, 'POINT_DATA %d\n', numel(scalar));
    fprintf(fid, 'SCALARS %s int 1\n', char(scalar_name));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%d\n', scalar(:));
    fclose(fid);
end

%% ==================== AGE FILTER ====================
function gw = parse_GA_from_filename(file_path)
    [~, name, ~] = fileparts(char(file_path));
    tok = regexp(name, '(?i)G[AW][_\-\.]?([0-9]+(\.[0-9]+)?)', 'tokens', 'once');
    if isempty(tok); gw = NaN; else; gw = str2double(tok{1}); end
end

function files_out = filter_files_by_GA(files_in, min_GA, max_GA)
    ga = arrayfun(@parse_GA_from_filename, files_in);
    bad = isnan(ga);
    if any(bad)
        fprintf('Warning: %d files have unreadable GA and were dropped.\n', nnz(bad));
    end
    keep = (~bad) & (ga >= min_GA) & (ga <= max_GA);
    files_out = files_in(keep);
    fprintf('GA filter: kept %d / %d files [%.2f, %.2f].\n', nnz(keep), numel(files_in), min_GA, max_GA);
end

%% ==================== LABELS / CLEANUP ====================
function labels = hard_labels_from_W(W)
    [~, labels] = max(W, [], 2);
    labels = int32(labels);
end

function labels_out = cleanup_labels_locked(labels_in, F, mask_cortex, cfg)
    labels_out = int32(labels_in(:));
    labels_out(~mask_cortex) = 0;

    do_component = any(strcmpi(cfg.island_cleanup_mode, ["component","both"]));
    do_majority  = any(strcmpi(cfg.island_cleanup_mode, ["majority","both"]));

    if do_component && cfg.min_island_size > 0
        labels_out = relabel_small_islands_locked(labels_out, F, mask_cortex, cfg.min_island_size);
    end
    if do_majority
        labels_out = relabel_by_neighbor_majority_locked(labels_out, F, mask_cortex, ...
            cfg.majority_vote_iters, cfg.majority_vote_frac);
    end
    labels_out(~mask_cortex) = 0;
end

function labels_out = relabel_small_islands_locked(labels_in, F, mask_cortex, minSize)
    labels_out = int32(labels_in(:));
    adj = build_vertex_adjacency(numel(labels_out), F);
    K = double(max(labels_out));

    for k = 1:K
        verts_k = find(labels_out == k & mask_cortex);
        visited = false(numel(labels_out), 1);
        visited(~mask_cortex | labels_out ~= k) = true;

        for v0 = verts_k(:)'
            if visited(v0); continue; end
            comp = bfs_component(v0, labels_out, k, adj, mask_cortex, visited);
            visited(comp) = true;

            if numel(comp) < minSize
                nb = unique(vertcat(adj{comp}));
                nb = nb(mask_cortex(nb));
                nb_lab = double(labels_out(nb));
                nb_lab = nb_lab(nb_lab ~= k & nb_lab ~= 0);
                if ~isempty(nb_lab)
                    cnt = accumarray(nb_lab, 1, [K,1], @sum, 0);
                    [~, newk] = max(cnt);
                    labels_out(comp) = int32(newk);
                end
            end
        end
    end
end

function comp = bfs_component(v0, labels, k, adj, mask_cortex, visited)
    q = v0;
    visited(v0) = true;
    comp = v0;
    while ~isempty(q)
        v = q(1); q(1) = [];
        nb = adj{v};
        nb = nb(mask_cortex(nb));
        nb = nb(labels(nb) == k & ~visited(nb));
        if ~isempty(nb)
            visited(nb) = true;
            q = [q; nb(:)]; %#ok<AGROW>
            comp = [comp; nb(:)]; %#ok<AGROW>
        end
    end
end

function labels_out = relabel_by_neighbor_majority_locked(labels_in, F, mask_cortex, iters, frac)
    labels_out = int32(labels_in(:));
    adj = build_vertex_adjacency(numel(labels_out), F);
    ctx = find(mask_cortex);
    K = double(max(labels_out));

    for t = 1:iters
        new_labels = labels_out;
        for i = ctx(:)'
            nb = adj{i};
            nb = nb(mask_cortex(nb));
            nb_lab = double(labels_out(nb));
            nb_lab = nb_lab(nb_lab ~= 0);
            if isempty(nb_lab); continue; end
            cnt = accumarray(nb_lab, 1, [K,1], @sum, 0);
            [mx, maj] = max(cnt);
            if mx / numel(nb_lab) >= frac
                new_labels(i) = int32(maj);
            end
        end
        labels_out = new_labels;
        labels_out(~mask_cortex) = 0;
    end
end

function adj = build_vertex_adjacency(nV, F)
    adj = cell(nV,1);
    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    for e = 1:size(E,1)
        i = E(e,1); j = E(e,2);
        adj{i}(end+1,1) = j;
        adj{j}(end+1,1) = i;
    end
end

%% ==================== STABILITY ====================
function instability = paper_style_instability(V_nmf, A_ctx, vertex_coords, k, cfg, labels_full, mesh, hemi, out_dir)
    nS = size(V_nmf,2);
    if nS < 2 * cfg.stab_min_scans
        instability = NaN;
        return;
    end

    d = nan(cfg.stab_nSplits, 1);
    for t = 1:cfg.stab_nSplits
        idx = randperm(nS);
        i1 = idx(1:floor(nS/2));
        i2 = idx(floor(nS/2)+1:end);

        [lab1, ~] = fit_partition_backend(V_nmf(:,i1), A_ctx, vertex_coords, k, cfg, "search");
        [lab2, ~] = fit_partition_backend(V_nmf(:,i2), A_ctx, vertex_coords, k, cfg, "search");

        if cfg.write_stability_vtk && should_write_stability_this_k(k, cfg)
            stab_dir = fullfile(out_dir, cfg.stability_vtk_dirname, sprintf("k%02d", k));
            if ~exist(stab_dir, 'dir'); mkdir(stab_dir); end

            [Vinf, Finf] = load_inflated_mesh(hemi, cfg, mesh);

            labels_full_1 = zeros(mesh.nV_full, 1, 'int32');
            labels_full_2 = zeros(mesh.nV_full, 1, 'int32');
            labels_full_1(mesh.mask_cortex) = lab1;
            labels_full_2(mesh.mask_cortex) = lab2;

            if cfg.cleanup_enable
                labels_full_1 = cleanup_labels_locked(labels_full_1, mesh.faces, mesh.mask_cortex, cfg);
                labels_full_2 = cleanup_labels_locked(labels_full_2, mesh.faces, mesh.mask_cortex, cfg);
            end

            write_vtk_point_scalar( ...
                fullfile(stab_dir, sprintf("labels_hemi_%s_k%02d_split%02d_half1.vtk", hemi, k, t)), ...
                Vinf, Finf, labels_full_1, "NMF_label");

            write_vtk_point_scalar( ...
                fullfile(stab_dir, sprintf("labels_hemi_%s_k%02d_split%02d_half2.vtk", hemi, k, t)), ...
                Vinf, Finf, labels_full_2, "NMF_label");
        end

        switch lower(cfg.stab_metric)
            case 'ari'
                d(t) = 0.5 * ((1 - adjusted_rand_index(lab1, labels_full)) + (1 - adjusted_rand_index(lab2, labels_full)));
            case 'nmi'
                d(t) = 0.5 * ((1 - normalized_mutual_information(lab1, labels_full)) + (1 - normalized_mutual_information(lab2, labels_full)));
            otherwise
                error('Unknown stability metric: %s', cfg.stab_metric);
        end
    end
    instability = mean(d, 'omitnan');
end

function score = split_half_stability(V_nmf, A_ctx, vertex_coords, k, cfg)
    nS = size(V_nmf,2);
    if nS < 2 * cfg.stab_min_scans
        score = NaN;
        return;
    end
    vals = nan(cfg.stab_nSplits,1);
    for t = 1:cfg.stab_nSplits
        idx = randperm(nS);
        i1 = idx(1:floor(nS/2));
        i2 = idx(floor(nS/2)+1:end);

        [lab1, ~] = fit_partition_backend(V_nmf(:,i1), A_ctx, vertex_coords, k, cfg, "search");
        [lab2, ~] = fit_partition_backend(V_nmf(:,i2), A_ctx, vertex_coords, k, cfg, "search");

        switch lower(cfg.stab_metric)
            case 'ari'; vals(t) = adjusted_rand_index(lab1, lab2);
            case 'nmi'; vals(t) = normalized_mutual_information(lab1, lab2);
            otherwise; error('Unknown stability metric: %s', cfg.stab_metric);
        end
    end
    score = mean(vals, 'omitnan');
end

function ari = adjusted_rand_index(L1, L2)
    L1 = grp2idx(categorical(double(L1(:))));
    L2 = grp2idx(categorical(double(L2(:))));
    n = numel(L1);
    C = accumarray([L1 L2], 1);
    nij = sum(C(:).*(C(:)-1)/2);
    ai = sum(C,2); bj = sum(C,1);
    a2 = sum(ai.*(ai-1)/2); b2 = sum(bj.*(bj-1)/2);
    expected = (a2*b2) / (n*(n-1)/2);
    maxval = 0.5*(a2+b2);
    den = maxval - expected;
    if abs(den) < 1e-12; ari = 0; else; ari = (nij - expected) / den; end
end

function nmi = normalized_mutual_information(L1, L2)
    L1 = grp2idx(categorical(double(L1(:))));
    L2 = grp2idx(categorical(double(L2(:))));
    C = accumarray([L1 L2], 1);
    Pij = C / sum(C(:));
    Pi = sum(Pij,2); Pj = sum(Pij,1);
    MI = 0;
    for i = 1:size(Pij,1)
        for j = 1:size(Pij,2)
            if Pij(i,j) > 0
                MI = MI + Pij(i,j) * log(Pij(i,j)/(Pi(i)*Pj(j)));
            end
        end
    end
    Hi = -sum(Pi(Pi>0).*log(Pi(Pi>0)));
    Hj = -sum(Pj(Pj>0).*log(Pj(Pj>0)));
    den = sqrt(Hi*Hj);
    if den < 1e-12; nmi = 0; else; nmi = MI / den; end
end

%% ==================== METRICS ====================
function err = compute_recon_error(V, R, type)
    switch lower(type)
        case 'abs_sse'
            err = sum(R(:).^2);
        case 'rel_fro'
            err = norm(R, 'fro') / (norm(V, 'fro') + eps);
        case 'rel_sse'
            err = (norm(R, 'fro')^2) / (norm(V, 'fro')^2 + eps);
        case 'col_rel_sse'
            err = mean(sum(R.^2, 1) ./ (sum(V.^2, 1) + eps));
        otherwise
            error('Unknown reconstruction error type: %s', type);
    end
end

function m = compute_mean_silhouette(X, labels, cfg)
    dist = 'correlation';
    if isfield(cfg, 'sil_distance') && ~isempty(cfg.sil_distance)
        dist = lower(cfg.sil_distance);
    end
    labels = labels(:);
    if numel(unique(labels)) < 2 || any(accumarray(labels,1) < 2)
        m = NaN;
        return;
    end
    try
        m = mean(silhouette(X, labels, dist), 'omitnan');
    catch
        m = NaN;
    end
end

function z = zscore_nan_safe(x)
    x = double(x(:));
    mu = mean(x, 'omitnan');
    sd = std(x, 'omitnan');
    if ~isfinite(sd) || sd < 1e-12
        z = zeros(size(x));
    else
        z = (x - mu) ./ sd;
    end
    z(~isfinite(z)) = 0;
end

%% ==================== OUTPUT ====================
function save_outputs(model, labels_ctx, labels_full, best_k, metrics, cfg, files, mesh, hemi, out_dir)
    W = [];
    H = [];
    if isfield(model, 'W'); W = model.W; end
    if isfield(model, 'H'); H = model.H; end

    save(fullfile(out_dir, sprintf('nmf_result_hemi_%s_k%02d.mat', hemi, best_k)), ...
        'W','H','labels_ctx','labels_full','best_k','metrics','cfg','files','mesh','-v7.3');

    writetable(table((1:numel(labels_full))', labels_full, ...
        'VariableNames', {'vertex_index','label'}), ...
        fullfile(out_dir, sprintf('labels_hemi_%s_k%02d.csv', hemi, best_k)));

    if cfg.write_vtk_labels
        [Vinf, Finf] = load_inflated_mesh(hemi, cfg, mesh);
        write_vtk_point_scalar(fullfile(out_dir, sprintf('labels_hemi_%s_k%02d.vtk', hemi, best_k)), ...
            Vinf, Finf, labels_full, 'NMF_label');
    end
end

function tf = should_write_sweep_vtk(cfg, mesh)
    tf = cfg.write_vtk_labels && isfield(mesh, 'faces') && ~isempty(mesh.faces) && ...
        any(strcmpi(cfg.vtk_mode, ["all","list"]));
end

function tf = should_write_this_k(k, cfg)
    switch lower(cfg.vtk_mode)
        case 'all'
            tf = true;
        case 'list'
            tf = ismember(k, cfg.vtk_k_list);
        otherwise
            tf = false;
    end
end

function shared_k = pick_shared_k(metricsL, metricsR, cfg)
    k = metricsL.k_range(:);
    joint = 0.5 * (metricsL.joint_score(:) + metricsR.joint_score(:));
    joint(k < cfg.k_min_select) = -Inf;
    [~, idx] = max(joint);
    shared_k = k(idx);
end

%% ==================== PLOTTING ====================
function plot_k_selection(metrics, hemi, out_dir)
    k = metrics.k_range;
    best_k = metrics.best_k;

    fig = figure('Color', 'w', 'Position', [100 100 1200 270]);
    tl = tiledlayout(1,4, 'Padding','compact', 'TileSpacing','compact'); %#ok<NASGU>

    plot_metric_panel(nexttile, k, metrics.silhouette_mean, best_k, ...
        'Silhouette coefficient', '(higher is better)');
    plot_metric_panel(nexttile, k, metrics.reconstruction_error, best_k, ...
        'Reconstruction error', '(lower is better)');
    plot_metric_panel(nexttile, k, metrics.instability, best_k, ...
        'Instability coefficient', '(lower is better)');
    plot_metric_panel(nexttile, k, metrics.joint_score, best_k, ...
        'Joint score', '(higher is better)');

    exportgraphics(fig, fullfile(out_dir, sprintf('k_selection_hemi_%s.png', hemi)), 'Resolution', 300);
    close(fig);
end

function plot_metric_panel(ax, k, y, best_k, ylab, subtitleText)
    axes(ax); hold(ax, 'on');
    plot(ax, k, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.1 0.1 0.1], 'Color', [0.1 0.1 0.1]);
    xline(ax, best_k, '--', 'LineWidth', 1.0, 'Color', [0.5 0.5 0.5]);

    idx = find(k == best_k, 1);
    if ~isempty(idx) && isfinite(y(idx))
        scatter(ax, best_k, y(idx), 45, 'filled', 'MarkerFaceColor', [0.1 0.1 0.1]);
        text(ax, best_k, y(idx), sprintf('  k = %d', best_k), ...
            'FontSize', 9, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end

    xlabel(ax, 'Region number k', 'FontName', 'Arial', 'FontSize', 10);
    ylabel(ax, ylab, 'FontName', 'Arial', 'FontSize', 10);
    title(ax, {ylab, subtitleText}, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');
    set(ax, 'Box', 'off', 'LineWidth', 0.8, 'FontName', 'Arial', 'FontSize', 9, 'TickDir', 'out');
    xlim(ax, [min(k) max(k)]);
end

function plot_k_selection_forced(out_dir, hemi, best_k)
    fig = figure('Color', 'w', 'Position', [100 100 420 120]);
    ax = axes(fig); axis(ax, 'off');
    text(ax, 0.05, 0.6, sprintf('Hemisphere %s\nForced k = %d', hemi, best_k), ...
        'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
    exportgraphics(fig, fullfile(out_dir, sprintf('forced_k_hemi_%s.png', hemi)), 'Resolution', 300);
    close(fig);
end

%% ==================== OPTIONAL SMOOTHING ====================
function Vsm = smooth_vertex_maps_masked_full(V, F, mask_cortex, iters)
    nV = size(V,1);
    adj = build_vertex_adjacency(nV, F);
    Vsm = V;
    for t = 1:iters
        Vnew = Vsm;
        for i = 1:nV
            nb = adj{i};
            if isempty(nb); continue; end
            if mask_cortex(i)
                nb = nb(mask_cortex(nb));
            else
                nb = nb(~mask_cortex(nb));
            end
            if isempty(nb); continue; end
            Vnew(i,:) = mean(Vsm([i; nb], :), 1);
        end
        Vsm = Vnew;
    end
end
%% write referene point location
function write_reference_points_vtk(out_path, vertices, faces, ref_idx)
    scalar = zeros(size(vertices,1),1);
    scalar(ref_idx) = 1;

    write_vtk_point_scalar(out_path, vertices, faces, scalar, "ref_points");
end

function [V_ctx, F_ctx] = extract_cortex_mesh(vertices, faces, mask)
    idx_map = zeros(size(mask));
    idx_map(mask) = 1:nnz(mask);

    keep_faces = mask(faces(:,1)) & mask(faces(:,2)) & mask(faces(:,3));
    F_ctx = faces(keep_faces, :);
    F_ctx = idx_map(F_ctx);

    V_ctx = vertices(mask, :);
end

function ref_idx = uniform_ref_indices_fps(vertices, nRef, mask, seed_idx)
    if nargin < 3 || isempty(mask)
        mask = true(size(vertices, 1), 1);
    end
    if nargin < 4
        seed_idx = [];
    end

    nV = size(vertices, 1);
    if size(vertices, 2) ~= 3
        error('vertices must be [nV x 3].');
    end
    if numel(mask) ~= nV
        error('mask must have length nV.');
    end

    valid_idx = find(mask);
    nValid = numel(valid_idx);

    if nValid == 0
        error('No valid vertices available for sampling.');
    end

    nRef = min(nRef, nValid);
    if nRef <= 0
        ref_idx = zeros(0,1);
        return;
    end

    X = double(vertices(valid_idx, :));

    % choose initial point
    if isempty(seed_idx)
        ctr = mean(X, 1);
        d0 = sum((X - ctr).^2, 2);
        [~, loc0] = max(d0);
    else
        loc0 = find(valid_idx == seed_idx, 1);
        if isempty(loc0)
            error('seed_idx is not inside the valid mask.');
        end
    end

    selected_local = zeros(nRef, 1);
    selected_local(1) = loc0;

    minDist2 = sum((X - X(loc0,:)).^2, 2);

    for t = 2:nRef
        [~, next_loc] = max(minDist2);
        selected_local(t) = next_loc;

        d2 = sum((X - X(next_loc,:)).^2, 2);
        minDist2 = min(minDist2, d2);
    end

    ref_idx = valid_idx(selected_local);
end