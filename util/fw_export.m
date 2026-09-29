function outfile = fw_export(fig, outfile, opts)
%FW_EXPORT  Write a figure to PNG (200 dpi), SVG or PDF (vector); creates the folder, prints the size.
%   outfile = fw_export(fig, outfile)
%   outfile = fw_export(fig, outfile, 'Resolution', 300, 'Verbose', false)
%   fig       figure handle, normally from fw_figure (its Paper* settings fix the pixel size;
%             the paper background is kept because fw_figure sets InvertHardcopy off)
%   outfile   target path; the extension selects the format: .png (raster at 'Resolution' dpi,
%             default fw_style().dpi = 200), .svg or .pdf (vector via print -vector)
%   'Verbose' true (default) prints '<file> written, <kB>'
%   Example:
%     fig = fw_figure(14, 9.8);  plot(1:10);  fw_export(fig, fullfile(tempdir, 'fw_test.png'))   % 1102 x 772 px

arguments
    fig (1,1) matlab.ui.Figure
    outfile (1,:) char
    opts.Resolution (1,1) double {mustBePositive} = 200
    opts.Verbose (1,1) logical = true
end

[folder, ~, ext] = fileparts(outfile);
if ~isempty(folder) && ~isfolder(folder)
    mkdir(folder);
end
switch lower(ext)
    case '.png'
        print(fig, outfile, '-dpng', sprintf('-r%d', round(opts.Resolution)));
    case '.svg'
        print(fig, outfile, '-dsvg', '-vector');
    case '.pdf'
        print(fig, outfile, '-dpdf', '-vector');
    otherwise
        error('fw_export:format', 'Unsupported extension "%s": use .png, .svg or .pdf', ext);
end
if opts.Verbose
    d = dir(outfile);
    fprintf('%s written, %.0f kB\n', outfile, d.bytes/1024);
end
end
