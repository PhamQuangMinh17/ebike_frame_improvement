clear; clc; close all;
inputFile = fullfile('ebike_topology_optimization_outputs','ebike_topology_optimization_step6.mat');
outDir = 'ebike_redesign_outputs';
cfg.thresholdDensity = 0.50;
cfg.autoTuneWidthToTopologyVolume = true;
cfg.minWidthScale = 0.75;
cfg.maxWidthScale = 2.50;
cfg.defaultWidthScaleIfNoTuning = 1.45;
cfg.baseWidth.chainStay_mm      = 45;
cfg.baseWidth.seatStay_mm       = 45;
cfg.baseWidth.seatTube_mm       = 50;
cfg.baseWidth.topTube_mm        = 58;
cfg.baseWidth.headTube_mm       = 60;
cfg.baseWidth.batteryRail_mm    = 42;
cfg.exportResolution = 300;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
if ~exist(inputFile, 'file')
    error(['Cannot find input file: %s\n' ...
           'Please run topology_optimization.m first, or check the folder name.'], inputFile);
end
S = load(inputFile);
requiredVars = {'feMesh','geometryData','baselineData','topOptData'};
for i = 1:numel(requiredVars)
    if ~isfield(S, requiredVars{i})
        error('Input MAT file does not contain required variable: %s', requiredVars{i});
    end
end
feMesh       = S.feMesh;
geometryData = S.geometryData;
baselineData = S.baselineData;
topOptData   = S.topOptData;
xCenters = feMesh.xCenters(:);
yCenters = feMesh.yCenters(:);
nelx = numel(xCenters);
nely = numel(yCenters);
[XC, YC] = ndgrid(xCenters, yCenters);
xPhysGrid = orientGrid(topOptData.xPhysGrid, nelx, nely, 'topOptData.xPhysGrid');
maskCode  = orientGrid(baselineData.maskCodeClean, nelx, nely, 'baselineData.maskCodeClean');
thresholdMask = xPhysGrid >= cfg.thresholdDensity;
if isfield(baselineData, 'elementVolume')
    elementVolume = baselineData.elementVolume;
else
    dx = mean(diff(xCenters));
    dy = mean(diff(yCenters));
    elementVolume = dx * dy * baselineData.material.thickness_mm;
end
if isfield(baselineData, 'elementArea')
    elementArea = baselineData.elementArea;
else
    dx = mean(diff(xCenters));
    dy = mean(diff(yCenters));
    elementArea = dx * dy;
end
baselineVolume = baselineData.materialVolumeBaseline;
targetTopologyVolume = topOptData.finalMaterialVolume;
P.RA     = rowPoint(geometryData.points.RA);
P.BB     = rowPoint(geometryData.points.BB);
P.ST     = rowPoint(geometryData.points.ST);
P.HT_top = rowPoint(geometryData.points.HT_top);
P.HT_bot = rowPoint(geometryData.points.HT_bot);
if isfield(geometryData.points, 'FA')
    P.FA = rowPoint(geometryData.points.FA);
else
    P.FA = [NaN NaN];
end
batteryPoly = ensureNx2(geometryData.batteryPoly, 'geometryData.batteryPoly');
batteryVoidMask = inpolygon(XC, YC, batteryPoly(:,1), batteryPoly(:,2));
batP1 = batteryPoly(1,:);
batP2 = batteryPoly(2,:);
batP3 = batteryPoly(3,:);
batP4 = batteryPoly(4,:);
memberDefs = struct('name', {}, 'path', {}, 'baseWidth_mm', {}, 'role', {}, 'startLabel', {}, 'endLabel', {});
memberDefs(end+1) = makeMember('M1 chain stay: RA to BB', ...
    [P.RA; P.BB], cfg.baseWidth.chainStay_mm, ...
    'Lower load path from rear support to motor/bottom-bracket mount.', 'RA', 'BB');
memberDefs(end+1) = makeMember('M2 seat stay: RA to ST', ...
    [P.RA; P.ST], cfg.baseWidth.seatStay_mm, ...
    'Diagonal load path from rider/seat region to rear support.', 'RA', 'ST');
memberDefs(end+1) = makeMember('M3 seat tube: BB to ST', ...
    [P.BB; P.ST], cfg.baseWidth.seatTube_mm, ...
    'Vertical/diagonal member connecting rider load to bottom-bracket motor region.', 'BB', 'ST');
memberDefs(end+1) = makeMember('M4 top tube: ST to HT top', ...
    [P.ST; P.HT_top], cfg.baseWidth.topTube_mm, ...
    'Upper load path toward the head tube, clearly retained by the density map.', 'ST', 'HT top');
memberDefs(end+1) = makeMember('M5 head tube: HT bottom to HT top', ...
    [P.HT_bot; P.HT_top], cfg.baseWidth.headTube_mm, ...
    'Front steering/head-tube support region; must remain structurally continuous.', 'HT bottom', 'HT top');
memberDefs(end+1) = makeMember('M6 upper battery rail', ...
    [P.BB; batP4; batP3; P.HT_bot], cfg.baseWidth.batteryRail_mm, ...
    'Upper rail around the battery cavity; keeps battery space open.', 'BB', 'HT bottom');
memberDefs(end+1) = makeMember('M7 lower battery rail', ...
    [P.BB; batP1; batP2; P.HT_bot], cfg.baseWidth.batteryRail_mm, ...
    'Lower rail around the battery cavity; keeps battery space open.', 'BB', 'HT bottom');
anchorDefs = struct('name', {}, 'point', {}, 'baseRadius_mm', {}, 'role', {});
anchorDefs(end+1) = makeAnchor('RA anchor', P.RA, getFieldOrDefault(geometryData.anchor, 'RA_radius', 45), ...
    'Rear axle/dropout support must not disappear.');
anchorDefs(end+1) = makeAnchor('BB / motor mount anchor', P.BB, getFieldOrDefault(geometryData.anchor, 'BB_radius', 60), ...
    'Bottom bracket and motor/pedalling load mount must remain solid.');
anchorDefs(end+1) = makeAnchor('ST anchor', P.ST, getFieldOrDefault(geometryData.anchor, 'ST_radius', 42), ...
    'Seat/rider load application region must remain solid.');
anchorDefs(end+1) = makeAnchor('HT top anchor', P.HT_top, getFieldOrDefault(geometryData.anchor, 'HT_radius', 45), ...
    'Head tube upper connection must remain solid.');
anchorDefs(end+1) = makeAnchor('HT bottom anchor', P.HT_bot, getFieldOrDefault(geometryData.anchor, 'HT_radius', 45), ...
    'Head tube lower connection and front support region must remain solid.');
if cfg.autoTuneWidthToTopologyVolume
    widthScale = tuneWidthScale(XC, YC, maskCode, batteryVoidMask, memberDefs, anchorDefs, ...
                                elementVolume, targetTopologyVolume, cfg.minWidthScale, cfg.maxWidthScale);
else
    widthScale = cfg.defaultWidthScaleIfNoTuning;
end
[redesignMask, memberMasks, anchorMasks] = buildRedesignMask(XC, YC, maskCode, batteryVoidMask, ...
    memberDefs, anchorDefs, widthScale);
redesignVolume = nnz(redesignMask) * elementVolume;
redesignArea   = nnz(redesignMask) * elementArea;
redesignRatioToTopology = redesignVolume / targetTopologyVolume;
redesignRatioToBaseline = redesignVolume / baselineVolume;
redesignReductionVsBaseline = 100 * (1 - redesignRatioToBaseline);
redesignMaskVector = vectorFromGridByElemIndex(redesignMask, feMesh.elemIndex, feMesh.nElements);
xPhysRedesignVector = double(redesignMaskVector);
memberName = cell(numel(memberDefs),1);
startLabel = cell(numel(memberDefs),1);
endLabel = cell(numel(memberDefs),1);
startX = zeros(numel(memberDefs),1);
startY = zeros(numel(memberDefs),1);
endX = zeros(numel(memberDefs),1);
endY = zeros(numel(memberDefs),1);
baseWidth_mm = zeros(numel(memberDefs),1);
finalWidth_mm = zeros(numel(memberDefs),1);
meanDensityInCorridor = zeros(numel(memberDefs),1);
role = cell(numel(memberDefs),1);
for i = 1:numel(memberDefs)
    pts = memberDefs(i).path;
    memberName{i} = memberDefs(i).name;
    startLabel{i} = memberDefs(i).startLabel;
    endLabel{i} = memberDefs(i).endLabel;
    startX(i) = pts(1,1);
    startY(i) = pts(1,2);
    endX(i) = pts(end,1);
    endY(i) = pts(end,2);
    baseWidth_mm(i) = memberDefs(i).baseWidth_mm;
    finalWidth_mm(i) = memberDefs(i).baseWidth_mm * widthScale;
    corridor = memberMasks(:,:,i) & (maskCode ~= 0);
    if any(corridor(:))
        meanDensityInCorridor(i) = mean(xPhysGrid(corridor));
    else
        meanDensityInCorridor(i) = NaN;
    end
    role{i} = memberDefs(i).role;
end
memberTable = table(memberName, startLabel, endLabel, startX, startY, endX, endY, ...
    baseWidth_mm, finalWidth_mm, meanDensityInCorridor, role);
writetable(memberTable, fullfile(outDir, 'redesign_members.csv'));
summaryItem = {'target_topology_volume_mm3'; ...
               'redesigned_binary_volume_mm3'; ...
               'baseline_full_solid_volume_mm3'; ...
               'redesign_volume_ratio_to_topology'; ...
               'redesign_volume_ratio_to_baseline'; ...
               'redesign_material_reduction_vs_baseline_percent'; ...
               'width_scale_used'; ...
               'number_of_redesigned_members'; ...
               'number_of_solid_elements_in_redesign'; ...
               'element_volume_mm3'};
summaryValue = [targetTopologyVolume; ...
                redesignVolume; ...
                baselineVolume; ...
                redesignRatioToTopology; ...
                redesignRatioToBaseline; ...
                redesignReductionVsBaseline; ...
                widthScale; ...
                numel(memberDefs); ...
                nnz(redesignMask); ...
                elementVolume];
summaryTable = table(summaryItem, summaryValue);
writetable(summaryTable, fullfile(outDir, 'redesign_summary.csv'));
fig = figure('Color','w', 'Position',[50 50 1600 900]);
imagesc(xCenters, yCenters, xPhysGrid');
set(gca, 'YDir','normal'); axis equal tight; hold on;
colormap(gca, flipud(gray(256))); cb = colorbar; ylabel(cb, 'Density x');
title('Topology density interpretation and retained load paths');
xlabel('x position (mm)'); ylabel('y position (mm)');
plotBattery(batteryPoly, 'Battery cavity / passive void');
plotAllMembers(memberDefs, widthScale, false);
plotKeyPointsAndBC(P, S.baselineData);
legend('Location','eastoutside');
exportFig(fig, fullfile(outDir, '21_topology_interpretation_load_paths.png'), cfg.exportResolution);
fig = figure('Color','w', 'Position',[50 50 1600 900]);
hold on; axis equal; grid on;
title('Final clean redesign: member centre-lines and design rules');
xlabel('x position (mm)'); ylabel('y position (mm)');
plotDesignDomainOutline(maskCode, xCenters, yCenters);
plotBattery(batteryPoly, 'Battery cavity');
plotAllMembers(memberDefs, widthScale, true);
plotAllAnchors(anchorDefs, widthScale);
plotKeyPointsAndBC(P, S.baselineData);
legend('Location','eastoutside');
exportFig(fig, fullfile(outDir, '22_redesigned_centerlines_and_rules.png'), cfg.exportResolution);
fig = figure('Color','w', 'Position',[50 50 1600 900]);
imagesc(xCenters, yCenters, double(redesignMask)');
set(gca, 'YDir','normal'); axis equal tight; hold on;
colormap(gca, gray(2)); cb = colorbar; ylabel(cb, '0 = void, 1 = solid');
title(sprintf('Redesigned binary geometry mask, volume ratio to baseline = %.3f', redesignRatioToBaseline));
xlabel('x position (mm)'); ylabel('y position (mm)');
plotBattery(batteryPoly, 'Battery cavity');
plotKeyPointsAndBC(P, S.baselineData);
exportFig(fig, fullfile(outDir, '23_redesigned_geometry_mask.png'), cfg.exportResolution);
fig = figure('Color','w', 'Position',[50 50 1600 900]);
imagesc(xCenters, yCenters, xPhysGrid');
set(gca, 'YDir','normal'); axis equal tight; hold on;
colormap(gca, flipud(gray(256))); cb = colorbar; ylabel(cb, 'Step 6 density x');
contour(xCenters, yCenters, double(redesignMask)', [0.5 0.5], 'r', 'LineWidth', 2.0);
contour(xCenters, yCenters, double(thresholdMask)', [0.5 0.5], 'b--', 'LineWidth', 1.0);
plotBattery(batteryPoly, 'Battery cavity');
title('Redesigned geometry overlay on topology density result');
xlabel('x position (mm)'); ylabel('y position (mm)');
legend({'redesigned geometry boundary', 'Step 6 threshold boundary'}, 'Location','eastoutside');
exportFig(fig, fullfile(outDir, '24_redesign_overlay_on_density.png'), cfg.exportResolution);
fig = figure('Color','w', 'Position',[50 50 1600 900]);
imagesc(xCenters, yCenters, 1 - double(redesignMask)');
set(gca, 'YDir','normal'); axis equal tight; hold on;
colormap(gca, gray(2));
plotBattery(batteryPoly, 'Battery cavity');
plotKeyPointsOnly(P);
title('Report-ready 2D redesigned e-bike frame geometry');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFig(fig, fullfile(outDir, '25_report_ready_redesigned_frame.png'), cfg.exportResolution);
fig = figure('Color','w', 'Position',[50 50 1600 900]);
comparisonGrid = zeros(nelx, nely);
comparisonGrid(thresholdMask & ~redesignMask) = 1;
comparisonGrid(redesignMask & ~thresholdMask) = 2;
comparisonGrid(redesignMask & thresholdMask)  = 3;
imagesc(xCenters, yCenters, comparisonGrid');
set(gca, 'YDir','normal'); axis equal tight; hold on;
colormap(gca, [1 1 1; 0.95 0.65 0.2; 0.4 0.7 1.0; 0.1 0.1 0.1]);
cb = colorbar;
set(cb, 'Ticks', [0 1 2 3], 'TickLabels', {'void', 'TO only', 'redesign only', 'overlap'});
plotBattery(batteryPoly, 'Battery cavity');
title('Comparison of optimized threshold topology and final clean redesign');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFig(fig, fullfile(outDir, '26_topology_vs_redesign_comparison.png'), cfg.exportResolution);
redesignData = struct();
redesignData.description = 'Final clean redesigned 2D e-bike frame geometry interpreted from the optimized density map';
redesignData.settings = cfg;
redesignData.widthScale = widthScale;
redesignData.memberDefs = memberDefs;
redesignData.anchorDefs = anchorDefs;
redesignData.memberTable = memberTable;
redesignData.summaryTable = summaryTable;
redesignData.redesignMask = redesignMask;
redesignData.redesignMaskVector = redesignMaskVector;
redesignData.xPhysRedesignVector = xPhysRedesignVector;
redesignData.thresholdMask = thresholdMask;
redesignData.batteryVoidMask = batteryVoidMask;
redesignData.batteryPoly = batteryPoly;
redesignData.redesignVolume_mm3 = redesignVolume;
redesignData.redesignArea_mm2 = redesignArea;
redesignData.targetTopologyVolume_mm3 = targetTopologyVolume;
redesignData.baselineVolume_mm3 = baselineVolume;
redesignData.redesignRatioToTopology = redesignRatioToTopology;
redesignData.redesignRatioToBaseline = redesignRatioToBaseline;
redesignData.redesignReductionVsBaseline_percent = redesignReductionVsBaseline;
save(fullfile(outDir, 'ebike_redesigned_geometry_step7.mat'), ...
    'redesignData', 'feMesh', 'geometryData', 'baselineData', 'topOptData', '-v7.3');
fprintf('\nFINAL CLEAN REDESIGN generated successfully.\n');
fprintf('Output folder: %s\n\n', fullfile(pwd, outDir));
fprintf('Interpretation summary:\n');
fprintf('  The Step 6 density map is converted into clean member-like geometry.\n');
fprintf('  Battery cavity is kept as passive void.\n');
fprintf('  RA, BB/motor mount, ST and HT regions are kept as passive solid anchors.\n\n');
fprintf('Redesign volume summary:\n');
fprintf('  Target Step 6 topology volume      = %.6e mm^3\n', targetTopologyVolume);
fprintf('  Redesigned binary geometry volume  = %.6e mm^3\n', redesignVolume);
fprintf('  Baseline full-solid volume         = %.6e mm^3\n', baselineVolume);
fprintf('  Redesign / topology volume ratio   = %.4f\n', redesignRatioToTopology);
fprintf('  Redesign / baseline volume ratio   = %.4f\n', redesignRatioToBaseline);
fprintf('  Material reduction vs baseline     = %.2f %%\n', redesignReductionVsBaseline);
fprintf('  Width scale used                   = %.4f\n\n', widthScale);
fprintf('Main redesigned members:\n');
for i = 1:numel(memberDefs)
    fprintf('  %-32s | width = %6.2f mm | mean Step 6 density in corridor = %.3f\n', ...
        memberDefs(i).name, memberDefs(i).baseWidth_mm * widthScale, meanDensityInCorridor(i));
end
fprintf('\nFiles exported:\n');
fprintf('  21_topology_interpretation_load_paths.png\n');
fprintf('  22_redesigned_centerlines_and_rules.png\n');
fprintf('  23_redesigned_geometry_mask.png\n');
fprintf('  24_redesign_overlay_on_density.png\n');
fprintf('  25_report_ready_redesigned_frame.png\n');
fprintf('  26_topology_vs_redesign_comparison.png\n');
fprintf('  redesign_members.csv\n');
fprintf('  redesign_summary.csv\n');
fprintf('  ebike_redesigned_geometry_step7.mat\n\n');
function M = orientGrid(M, nelx, nely, name)
    if isequal(size(M), [nelx, nely])
        return;
    elseif isequal(size(M), [nely, nelx])
        M = M';
    else
        error('%s has size %dx%d, expected %dx%d.', name, size(M,1), size(M,2), nelx, nely);
    end
end
function p = rowPoint(p)
    p = p(:).';
    if numel(p) ~= 2
        error('Point must contain exactly two coordinates.');
    end
end
function P = ensureNx2(P, name)
    if size(P,2) == 2
        return;
    elseif size(P,1) == 2
        P = P';
    else
        error('%s must be Nx2 or 2xN.', name);
    end
end
function val = getFieldOrDefault(S, fieldName, defaultValue)
    if isfield(S, fieldName)
        val = S.(fieldName);
    else
        val = defaultValue;
    end
end
function m = makeMember(name, path, baseWidth_mm, role, startLabel, endLabel)
    m.name = name;
    m.path = path;
    m.baseWidth_mm = baseWidth_mm;
    m.role = role;
    m.startLabel = startLabel;
    m.endLabel = endLabel;
end
function a = makeAnchor(name, point, baseRadius_mm, role)
    a.name = name;
    a.point = point;
    a.baseRadius_mm = baseRadius_mm;
    a.role = role;
end
function widthScale = tuneWidthScale(XC, YC, maskCode, batteryVoidMask, memberDefs, anchorDefs, ...
                                     elementVolume, targetVolume, minScale, maxScale)
    low = minScale;
    high = maxScale;
    bestScale = low;
    bestErr = inf;
    for iter = 1:40
        mid = 0.5 * (low + high);
        mask = buildRedesignMask(XC, YC, maskCode, batteryVoidMask, memberDefs, anchorDefs, mid);
        vol = nnz(mask) * elementVolume;
        err = abs(vol - targetVolume);
        if err < bestErr
            bestErr = err;
            bestScale = mid;
        end
        if vol < targetVolume
            low = mid;
        else
            high = mid;
        end
    end
    widthScale = bestScale;
end
function [solidMask, memberMasks, anchorMasks] = buildRedesignMask(XC, YC, maskCode, batteryVoidMask, ...
                                                                   memberDefs, anchorDefs, widthScale)
    solidMask = false(size(XC));
    memberMasks = false([size(XC), numel(memberDefs)]);
    anchorMasks = false([size(XC), numel(anchorDefs)]);
    for i = 1:numel(memberDefs)
        width = memberDefs(i).baseWidth_mm * widthScale;
        memberMasks(:,:,i) = maskFromPolyline(XC, YC, memberDefs(i).path, width);
        solidMask = solidMask | memberMasks(:,:,i);
    end
    for i = 1:numel(anchorDefs)
        radius = anchorDefs(i).baseRadius_mm * widthScale;
        p = anchorDefs(i).point;
        anchorMasks(:,:,i) = sqrt((XC - p(1)).^2 + (YC - p(2)).^2) <= radius;
        solidMask = solidMask | anchorMasks(:,:,i);
    end
    solidMask = solidMask & (maskCode ~= 0);
    solidMask(batteryVoidMask) = false;
    solidMask(maskCode == 2) = true;
end
function m = maskFromPolyline(XC, YC, path, width)
    m = false(size(XC));
    for k = 1:(size(path,1)-1)
        d = distanceToSegment(XC, YC, path(k,:), path(k+1,:));
        m = m | (d <= width/2);
    end
end
function d = distanceToSegment(X, Y, p1, p2)
    vx = p2(1) - p1(1);
    vy = p2(2) - p1(2);
    L2 = vx^2 + vy^2;
    if L2 < eps
        d = sqrt((X - p1(1)).^2 + (Y - p1(2)).^2);
        return;
    end
    t = ((X - p1(1)) .* vx + (Y - p1(2)) .* vy) ./ L2;
    t = max(0, min(1, t));
    projX = p1(1) + t .* vx;
    projY = p1(2) + t .* vy;
    d = sqrt((X - projX).^2 + (Y - projY).^2);
end
function vec = vectorFromGridByElemIndex(maskGrid, elemIndex, nElements)
    elemIndex = round(elemIndex);
    if ~isequal(size(maskGrid), size(elemIndex))
        if isequal(size(maskGrid), fliplr(size(elemIndex)))
            maskGrid = maskGrid';
        else
            error('maskGrid and elemIndex have incompatible sizes.');
        end
    end
    vec = false(1, nElements);
    for i = 1:size(elemIndex,1)
        for j = 1:size(elemIndex,2)
            eID = elemIndex(i,j);
            if eID >= 1 && eID <= nElements
                vec(eID) = maskGrid(i,j);
            end
        end
    end
end
function plotBattery(batteryPoly, displayName)
    patch('XData', batteryPoly(:,1), 'YData', batteryPoly(:,2), ...
          'FaceColor', [1 0.2 0.2], 'FaceAlpha', 0.12, ...
          'EdgeColor', [1 0 0], 'LineWidth', 2.0, 'DisplayName', displayName);
end
function plotAllMembers(memberDefs, widthScale, showWidth)
    for i = 1:numel(memberDefs)
        pts = memberDefs(i).path;
        if showWidth
            lw = max(2, memberDefs(i).baseWidth_mm * widthScale / 8);
        else
            lw = 2.2;
        end
        plot(pts(:,1), pts(:,2), '-', 'LineWidth', lw, 'DisplayName', memberDefs(i).name);
        midID = max(1, round(size(pts,1)/2));
        text(pts(midID,1), pts(midID,2), sprintf(' M%d', i), ...
             'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);
    end
end
function plotAllAnchors(anchorDefs, widthScale)
    th = linspace(0, 2*pi, 120);
    for i = 1:numel(anchorDefs)
        p = anchorDefs(i).point;
        r = anchorDefs(i).baseRadius_mm * widthScale;
        plot(p(1) + r*cos(th), p(2) + r*sin(th), 'k-', 'LineWidth', 1.5, ...
             'DisplayName', anchorDefs(i).name);
    end
end
function plotKeyPointsAndBC(P, baselineData)
    names = {'RA','BB','ST','HT top','HT bottom'};
    pts = [P.RA; P.BB; P.ST; P.HT_top; P.HT_bot];
    plot(pts(:,1), pts(:,2), 'ko', 'MarkerFaceColor','y', 'MarkerSize',7, 'DisplayName','key points');
    for i = 1:size(pts,1)
        text(pts(i,1)+10, pts(i,2)+8, names{i}, 'FontWeight','bold', 'Color','k', 'BackgroundColor','w');
    end
    if isfield(baselineData, 'loadSpecs')
        try
            for i = 1:numel(baselineData.loadSpecs)
                p = rowPoint(baselineData.loadSpecs(i).targetPoint);
                Fx = baselineData.loadSpecs(i).Fx;
                Fy = baselineData.loadSpecs(i).Fy;
                scale = 0.10;
                quiver(p(1), p(2), Fx*scale, Fy*scale, 0, 'r', 'LineWidth',2, 'MaxHeadSize',2, ...
                       'DisplayName','loads');
            end
        catch
        end
    end
end
function plotKeyPointsOnly(P)
    names = {'RA','BB','ST','HT top','HT bottom'};
    pts = [P.RA; P.BB; P.ST; P.HT_top; P.HT_bot];
    plot(pts(:,1), pts(:,2), 'ro', 'MarkerFaceColor','r', 'MarkerSize',6);
    for i = 1:size(pts,1)
        text(pts(i,1)+10, pts(i,2)+8, names{i}, 'FontWeight','bold', 'Color','r', 'BackgroundColor','w');
    end
end
function plotDesignDomainOutline(maskCode, xCenters, yCenters)
    contour(xCenters, yCenters, double(maskCode' ~= 0), [0.5 0.5], 'k--', 'LineWidth', 1.2, ...
            'DisplayName','original design domain');
end
function exportFig(figHandle, outPath, resolution)
    try
        exportgraphics(figHandle, outPath, 'Resolution', resolution);
    catch
        saveas(figHandle, outPath);
    end
end
