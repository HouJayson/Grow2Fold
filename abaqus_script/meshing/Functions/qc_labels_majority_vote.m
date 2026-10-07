function labelsOut = qc_labels_majority_vote(elements, labelsIn, varargin)
%QC_LABELS_MAJORITY_VOTE_KEEP_MAIN_ZERO
% Keeps the largest connected component of every label, including label 0.
% Smaller disconnected components are reassigned by neighbor majority vote.
%
% This is useful when label 0 is a real large region, but isolated 0 islands
% are misassigned regions.

% ---- options ----
p = inputParser;
p.addParameter('badLabels', -10);
p.addParameter('zeroLabel', 0);
p.addParameter('minSharedNodes', 1);
p.addParameter('maxIters', 5);          % repeat correction several times
p.addParameter('minIslandSize', 0);     % if >0, only fix components <= this size
p.parse(varargin{:});

badLabels      = p.Results.badLabels;
zeroLabel      = p.Results.zeroLabel;
minSharedNodes = p.Results.minSharedNodes;
maxIters       = p.Results.maxIters;
minIslandSize  = p.Results.minIslandSize;

labelsOut = labelsIn(:);

nElem = size(elements, 1);
nen   = size(elements, 2);
nNode = max(elements(:));

% ---- build element adjacency ----
rowsAll = double(elements(:));
colsAll = double(repmat((1:nElem)', nen, 1));

IncAll = sparse(rowsAll, colsAll, 1, double(nNode), double(nElem));

Aall = IncAll' * IncAll;
Aall = Aall - spdiags(diag(Aall), 0, nElem, nElem);

AdjAll = Aall >= minSharedNodes;

% ---- majority-vote correction ----
for it = 1:maxIters

    changed = false;

    labelsToFix = setdiff(unique(labelsOut)', badLabels);

    for L = labelsToFix

        elemIdxL = find(labelsOut == L);
        if numel(elemIdxL) <= 1
            continue;
        end

        AdjSub = AdjAll(elemIdxL, elemIdxL);
        Gsub   = graph(AdjSub);

        compIdx = conncomp(Gsub);
        nComp   = max(compIdx);

        if nComp <= 1
            continue;
        end

        compCounts = accumarray(compIdx(:), 1);

        % keep largest connected component
        [~, largestCompID] = max(compCounts);

        for c = 1:nComp

            if c == largestCompID
                continue;
            end

            compSize = compCounts(c);

            % optional: only correct small islands
            if minIslandSize > 0 && compSize > minIslandSize
                continue;
            end

            compGlobalIDs = elemIdxL(compIdx == c);

            neighMask = any(AdjAll(compGlobalIDs, :), 1);
            neighIDs  = find(neighMask);
            neighIDs  = setdiff(neighIDs, compGlobalIDs);

            neighLabels = labelsOut(neighIDs);

            % exclude bad labels and the current label itself
            neighLabels = neighLabels(~ismember(neighLabels, [badLabels, L]));

            if isempty(neighLabels)
                continue;
            end

            newLabel = mode_with_smallest_tie(neighLabels);

            labelsOut(compGlobalIDs) = newLabel;
            changed = true;
        end
    end

    if ~changed
        break;
    end
end

end

% ----- helper -----
function m = mode_with_smallest_tie(x)
u = unique(x);
cnt = zeros(numel(u), 1);

for k = 1:numel(u)
    cnt(k) = sum(x == u(k));
end

mx = max(cnt);
m = min(u(cnt == mx));
end