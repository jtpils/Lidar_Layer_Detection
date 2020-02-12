function [layerInfo, PD, PN] = VDE_cld(signal, height, BG, minLayerDepth, ...
                                       minHeight)
%VDE_CLD cloud layer detection with VDE method. THis method only required elstic
%signal.
%
%Inputs:
%   signal: array
%       raw signal without background. [photon count]
%   height: array
%       height above the ground. [km]
%   BG: numeric
%       background signal. [photon count]
%   minLayerDepth: numeric
%       minimun layer geometrical depth (default: 0.2). [km]
%   minHeight: numeric
%       minimum height to start the searching with (default: 0.4). [km]
%
%Returns:
%   layerInfo: struct array
%       id: numeric
%           identity of the layer.
%       baseHeight: numeric
%           the layer base height. [km]
%       topHeight: numeric
%           the layer top height. [km]
%       layerDepth: numeric
%           geometrical depth of the layer. [km]
%       flagCloud: logical
%           cloud flag.
%   PD: array
%       SDP signal.
%   PN: array
%       VDE signal
%
%References:
%   1. Zhao, C., Y. Wang, Q. Wang, Z. Li, Z. Wang, and D. Liu (2014), A new
%      cloud and aerosol layer detection method based on micropulse lidar 
%      measurements, Journal of Geophysical Research: Atmospheres, 119(11),
%      6788-6802.
%
%History:
%   2020-02-11. First edition by Zhenping.

if nargin < 3
    error('Not enough inputs.');
end

if ~ exist('minHeight', 'var')
    minHeight = 0.4;
end

if ~ exist('minLayerDepth', 'var')
    minLayerDepth = 0.2;
end

layerInfo = struct('id', {}, 'baseHeight', {}, 'topHeight', {}, ...
                   'layerDepth', {}, 'flagCloud', {});

minIndex = find(height >= minHeight, 1);
if isempty(minIndex)
    warning('minHeight is out of range.');
    return;
end

P = signal(minIndex:end);
height = height(minIndex:end);

%% 1. Semi-Discretization Process (SDP)
% noise_level = sqrt(BG) * 3;
noise_level = sqrt(BG + P);
Ps = transpose(smooth(P, max([ceil(0.09/(height(2) - height(1))), 3]), 'moving'));

% bottom to top semi-discretization
PD1 = Ps;
for indx = 2:length(PD1)
    if abs(PD1(indx) - PD1(indx - 1)) < max([noise_level(indx) * 3, sqrt(BG) * 3])
        PD1(indx) = PD1(indx - 1);
    end
end

% top to bottom semi-discretization
PD2 = Ps;
for indx = (length(PD2) - 1):-1:1
    if abs(PD2(indx) - PD2(indx + 1)) < max([noise_level(indx) * 3, sqrt(BG) * 3])
        PD2(indx) = PD2(indx + 1);
    end
end

PD = mean([PD1; PD2], 1);

%% 2. Value Distribution Equalization (VDE) Process
[RS, IS] = sort(PD);
MA = RS(end);
MI = RS(1);
PE = (1:length(RS)) / length(RS);
epsilon = 1e-6;

for indx = 2:length(RS)
    if abs(RS(indx) - RS(indx - 1)) <= epsilon
        PE(indx) = PE(indx - 1);
    end
end

yi = PE .* (MA - MI) + MI;

[~, RIS] = sort(IS);
PN = yi(RIS);

%% 3. Layer detection
BZ = (length(PN):-1:1)/length(PN) * (MA - MI) + MI;
layerN = 0;
[L, nLayer] = bwlabel(PN > BZ);
for iLayer = 1:nLayer
    baseIndex = find(L == iLayer, 1);
    topIndex = find(L == iLayer, 1, 'last');

    layerDepth = height(topIndex) - height(baseIndex);

    if (layerDepth >= minLayerDepth) && ...
       (mean(P(L == iLayer)) >= max([mean(noise_level(L == iLayer)*2), sqrt(BG)*3]))
        layerN = layerN + 1;
        layerInfo(layerN).id = layerN;
        layerInfo(layerN).baseHeight = height(baseIndex);
        layerInfo(layerN).topHeight = height(topIndex);
        layerInfo(layerN).layerDepth = layerDepth;
        layerInfo(layerN).flagCloud = false;
    end
end

%% 4. Layer classification
for iLayer = 1:length(layerInfo)

    indx = (height >= layerInfo(iLayer).baseHeight) & ...
           (height <= layerInfo(iLayer).topHeight);
    sig = Ps(indx) .* height(indx).^2;
    sig(sig<=0) = NaN;
    Fz = diff(log(sig .* height(indx).^2)) ./ diff(height(indx));

    T = max(sig) / sig(1);
    D = min(Fz);
    layerHeight = mean([height(indx(1)), height(indx(end))]);

    if layerHeight <= 5
        if T > 4 || D < -7
            layerInfo(iLayer).flagCloud = true;
        end
    else
        if T > 1.5 || D < -7
            layerInfo(iLayer).flagCloud = true;
        end
    end
end

end