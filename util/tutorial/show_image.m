function fig = show_image(file, opts)
%SHOW_IMAGE  Show a PNG from the MATLAB path in a figure of the standard width (TUTORIAL.m section (4)).
%   fig = show_image('DLM_FEM_Coupling.png')
%   file     file name of an image on the path (startup.m puts docs/figures there), or a full path
%   The figure is fw_style().size.standard(1) wide, its height follows the aspect ratio of the
%   image, and the image fills it without axes.
%   Options (defaults)
%     'Name'   file   figure name
%   fig      the figure
%   See also fw_figure.

arguments
    file (1,:) char
    opts.Name (1,:) char = file
end

S = fw_style();
path_img = which(file);
if isempty(path_img), path_img = file; end
img = imread(path_img);
fig = fw_figure(S.size.standard(1), S.size.standard(1)*size(img,1)/size(img,2), 'Name', opts.Name);
ax  = axes(fig, 'Position', [0 0 1 1]);
image(ax, img); axis(ax, 'image', 'off')
end
