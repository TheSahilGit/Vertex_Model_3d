function [files_out, its_out] = select_snapshots(files, its, it_start, it_end, it_stride)
% SELECT_SNAPSHOTS  Filter/subsample a (files, its) list (as returned
% by list_snapshots.m) down to an inclusive iteration range and a
% stride, all expressed in terms of the SAVED snapshot iteration
% numbers (i.e. the it_dumps cadence from para_Simulation.dat), not
% raw simulation steps and not a plain index into the file list.
%
%   [f,i] = select_snapshots(files, its, 0, inf, 1)     % every saved snapshot
%   [f,i] = select_snapshots(files, its, 1000, 4000, 2) % every OTHER saved
%                                                        % snapshot with
%                                                        % 1000 <= it <= 4000
%
% it_start/it_end may be -inf/inf to leave that end of the range open.
% it_stride subsamples the (already range-filtered) list, e.g.
% it_stride=1 keeps every entry, it_stride=3 keeps every third one.

if nargin < 3 || isempty(it_start),  it_start  = -inf; end
if nargin < 4 || isempty(it_end),    it_end    = inf;  end
if nargin < 5 || isempty(it_stride), it_stride = 1;    end
if it_stride < 1
    error('select_snapshots:stride', 'it_stride must be >= 1');
end

mask = (its >= it_start) & (its <= it_end);
files = files(mask);
its   = its(mask);

[its_sorted, order] = sort(its);
files = files(order);

idx = 1:round(it_stride):numel(files);
files_out = files(idx);
its_out   = its_sorted(idx);
end
