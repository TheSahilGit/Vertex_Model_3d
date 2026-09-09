function [files_out, its_out] = select_snapshots(files, its, wanted_its)
% SELECT_SNAPSHOTS  Pick out the snapshots at the requested SAVED
% iteration numbers (as returned by list_snapshots.m), e.g.:
%
%   select_snapshots(files, its, 1000:100:5000)   % that whole range
%   select_snapshots(files, its, 5000)            % a single iteration
%                                                  % (a scalar works
%                                                  % exactly like a
%                                                  % 1-element array)
%   select_snapshots(files, its)                  % every available
%                                                  % snapshot (wanted_its
%                                                  % omitted/empty)
%
% Each requested iteration is matched to the nearest AVAILABLE one
% (only some iterations were actually saved, per it_dumps); duplicate
% matches are dropped, and the result is always sorted into
% chronological order regardless of the input order, since that's the
% order a video should play in.

if nargin < 3 || isempty(wanted_its)
    [its_out, order] = sort(its);
    files_out = files(order);
    return;
end

wanted_its = wanted_its(:)';  % row vector, whatever shape/orientation came in

idx = zeros(size(wanted_its));
for k = 1:numel(wanted_its)
    [d, idx(k)] = min(abs(its - wanted_its(k)));
    if d > 0
        fprintf('select_snapshots: requested it=%g not saved, using nearest available it=%g instead.\n', ...
                wanted_its(k), its(idx(k)));
    end
end
idx = unique(idx, 'stable');

its_sel = its(idx);
[its_out, order] = sort(its_sel);
files_out = files(idx(order));
end
