function [files, its] = list_diagnostics(data_dir)
% LIST_DIAGNOSTICS  List available diag_<it>.dat files in chronological
% order, preferring data/diag_list.txt (written incrementally by the
% Fortran code) and falling back to a directory glob if that manifest
% is missing. Also returns the saved iteration number of each file
% (parsed from its 'diag_<it>.dat' filename). Mirrors
% list_snapshots.m exactly, for the diagnostics series instead of the
% vertex snapshots.
%
%   files = list_diagnostics()                 % looks in 'data/'
%   [files, its] = list_diagnostics('data')

if nargin < 1
    data_dir = 'data';
end

manifest = fullfile(data_dir, 'diag_list.txt');
files = {};
if exist(manifest, 'file')
    fid = fopen(manifest, 'r');
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);
        if isempty(line) || line(1) == '%'
            continue;
        end
        parts = strsplit(line);
        if numel(parts) >= 2
            % Re-root just the basename onto the caller's data_dir
            % (see list_snapshots.m for why: the manifest stores the
            % path the Fortran executable itself used).
            [~, base, ext] = fileparts(parts{2});
            files{end+1} = fullfile(data_dir, [base ext]); %#ok<AGROW>
        end
    end
    fclose(fid);
end

if isempty(files)
    d = dir(fullfile(data_dir, 'diag_*.dat'));
    [~, order] = sort({d.name});
    d = d(order);
    files = fullfile({d.folder}, {d.name});
end

its = zeros(1, numel(files));
for k = 1:numel(files)
    [~, base] = fileparts(files{k});
    tok = regexp(base, 'diag_(\d+)', 'tokens', 'once');
    if isempty(tok)
        its(k) = NaN;
    else
        its(k) = str2double(tok{1});
    end
end
end
