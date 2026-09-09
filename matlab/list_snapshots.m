function files = list_snapshots(data_dir)
% LIST_SNAPSHOTS  List available snapshot files in chronological order,
% preferring data/dump_list.txt (written incrementally by the Fortran
% code) and falling back to a directory glob if that manifest is
% missing.
%
%   files = list_snapshots()            % looks in 'data/'
%   files = list_snapshots('data')

if nargin < 1
    data_dir = 'data';
end

manifest = fullfile(data_dir, 'dump_list.txt');
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
            % The manifest stores the path the Fortran executable itself
            % used (typically 'data/snap_...dat', relative to ITS OWN
            % working directory at run time). Re-root just the basename
            % onto the caller's data_dir instead of trusting that
            % embedded path, so this works regardless of the caller's
            % current working directory.
            [~, base, ext] = fileparts(parts{2});
            files{end+1} = fullfile(data_dir, [base ext]); %#ok<AGROW>
        end
    end
    fclose(fid);
end

if isempty(files)
    d = dir(fullfile(data_dir, 'snap_*.dat'));
    [~, order] = sort({d.name});
    d = d(order);
    files = fullfile({d.folder}, {d.name});
end
end
