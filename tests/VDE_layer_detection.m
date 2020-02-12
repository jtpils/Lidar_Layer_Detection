clc; close all;

%% Add searching path
projectDir = fileparts(fileparts(mfilename('fullpath')));
libDir = fullfile(projectDir, 'lib');
addpath(libDir)
addincludepath;

%% initialization
% meteorFile = fullfile(projectDir, 'data', 'radiosonde_57494_20140101_0000.nc');
tRange = [datenum(2013, 8, 1, 0, 40, 0), datenum(2013, 8, 1, 23, 45, 0)];
hRange = [0, 20];
dataFile = fullfile(projectDir, 'data', sprintf('MUA_PLidar_%s.h5', datestr(tRange(1), 'yyyymmdd')));
figFile = fullfile(projectDir, 'img', sprintf('VDE_layer_detection.png'));

%% read data
PCR = h5read(dataFile, '/CH1/DataProcessed/LicelGluedData_PC');
time = floor(tRange(1)) + (1:size(PCR, 1)) * datenum(0, 1, 0, 0, 1, 0);
height = (1:size(PCR, 2)) * 0.03 + 0.074;
flag = (time <= tRange(2)) & (time >= tRange(1));
PCR = PCR(flag, :);
time = time(flag);
minHeight = 1.5;   % minimum height for searching layers. Unit: km
minLayerDepth = 0.1;   % minimun geometrical depth of the layer. Unit: km.

%% searching the layers
layer_status = zeros(size(PCR));   % 0: unknown; 1: aerosol; 2: cloud
for iTime = 1:sum(flag)
    PC = PCR(iTime, :) * 200;
    BG = mean(PC(1500:2000));
    PC = PC - BG;

    [layerInfo, PD, PN] = VDE_cld(PC, height, BG, minLayerDepth, minHeight);

    for iLayer = 1:length(layerInfo)
        layer_index = (height >= layerInfo(iLayer).baseHeight) & ...
                      (height <= layerInfo(iLayer).topHeight);
        if layerInfo(iLayer).flagCloud
            layer_status(iTime, layer_index) = 2;
        else
            layer_status(iTime, layer_index) = 1;
        end
    end
end

%% data visualization
load('myjet_colormap.mat');

figure('Position', [0, 20, 600, 500], 'Units', 'Pixels');
figPos = subfigPos([0.1, 0.13, 0.73, 0.8], 2, 1, 0, 0.03);

subplot('Position', figPos(1, :), 'Units', 'Normalized');

RCS = (PCR - repmat(mean(PCR(:, 1500:2000), 2), 1, size(PCR, 2))) .* (repmat(height, size(PCR, 1), 1)).^2;
p1 = pcolor(time, height, transpose(RCS)); hold on;
p1.EdgeColor = 'None';

ylabel('Height (km)');

set(gca, 'XTickLabel', '', 'XMinorTick', 'on', ...
         'YTick', linspace(hRange(1), hRange(2), 6), 'YMinorTick', 'on', ...
         'Box', 'on', 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 12);

text(0.02, 0.8, '(a) Range-corrected signal at 532 nm', 'Units', 'Normalized', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'Bold');

xlim(tRange);
ylim(hRange);
caxis([0, 1e3]);
colormap(gca, myjet);

cb = colorbar('Position', [figPos(1, 1) + figPos(1, 3) + 0.03, figPos(1, 2) + 0.05, 0.02, figPos(1, 4) - 0.10], 'Units', 'Normalized');
set(gca, 'TickDir', 'out', 'Box', 'on');
titleHandle = get(cb, 'Title');
set(titleHandle, 'string', '[a.u.]', 'FontSize', 12);

subplot('Position', figPos(2, :), 'Units', 'Normalized');
load('layer_status_colormap.mat');
p1 = pcolor(time, height, transpose(layer_status)); hold on;

p1.EdgeColor = 'None';

xlabel('Time (UTC)');
ylabel('Height (km)');

set(gca, 'XMinorTick', 'on', ...
         'YTick', linspace(hRange(1), hRange(2), 6), 'YMinorTick', 'on', ...
         'Box', 'on', 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 12);
datetick(gca, 'x', 'HH:MM', 'keepticks');


text(0.02, 0.8, '(b) layer status', 'Units', 'Normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'Bold');

xlim(tRange);
ylim(hRange);
caxis([0, 2]);

tickLabels = {'Unknown', ...
              'aerosol', ...
              'cloud'};
cb = colorbar('Position', [figPos(2, 1) + figPos(2, 3) + 0.03, figPos(2, 2) + 0.05, 0.02, figPos(2, 4) - 0.10], 'Units', 'Normalized');
colormap(gca, layer_status_colormap);
set(cb, 'TickDir', 'out', 'Box', 'on');
set(cb, 'ytick', (0.5:1:2.5)/3*2, 'yticklabel', tickLabels);

set(findall(gcf, '-Property', 'FontName'), 'FontName', 'Times New Roman');

fprintf('Export figure to %s\n', figFile);
export_fig(gcf, figFile, '-r300');
