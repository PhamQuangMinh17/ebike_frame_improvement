clear; clc; close all;
meshPreset = 'debug';
customNelx = 120;
customNely = 60;
loads.rider_Fx_N   = 0;
loads.rider_Fy_N   = -800;
loads.battery_Fx_N = 0;
loads.battery_Fy_N = -100;
loads.motor_Fx_N   = 120;
loads.motor_Fy_N   = -150;
plotFullGridLines = true;
maxGridLinesToPlot = 350;
geometryFile = fullfile(pwd, 'ebike_geometry_outputs', 'ebike_geometry_yasin2023.mat');
outDir = fullfile(pwd, 'ebike_fe_mesh_outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
if ~exist(geometryFile, 'file')
    error(['Cannot find geometry file:\n  %s\n\n', ...
           'Run build_ebike_geometry.m first, and make sure the folder ', ...
           'ebike_geometry_outputs is in the current MATLAB folder.'], geometryFile);
end
load(geometryFile, 'geometryData');
P = geometryData.points;
switch lower(meshPreset)
    case 'debug'
        nelx = 120;
        nely = 60;
    case 'medium'
        nelx = 180;
        nely = 90;
    case 'final'
        nelx = 220;
        nely = 110;
    case 'custom'
        nelx = customNelx;
        nely = customNely;
    otherwise
        error('Unknown meshPreset. Use debug, medium, final, or custom.');
end
xmin = geometryData.grid.xEdges(1);
xmax = geometryData.grid.xEdges(end);
ymin = geometryData.grid.yEdges(1);
ymax = geometryData.grid.yEdges(end);
xEdges = linspace(xmin, xmax, nelx + 1);
yEdges = linspace(ymin, ymax, nely + 1);
xCenters = 0.5 * (xEdges(1:end-1) + xEdges(2:end));
yCenters = 0.5 * (yEdges(1:end-1) + yEdges(2:end));
nodeID = reshape(1:(nelx + 1) * (nely + 1), nely + 1, nelx + 1);
nodeX = repmat(xEdges, nely + 1, 1);
nodeY = repmat(yEdges(:), 1, nelx + 1);
nNodes = (nelx + 1) * (nely + 1);
nodeCoords = zeros(nNodes, 2);
nodeCoords(nodeID(:), 1) = nodeX(:);
nodeCoords(nodeID(:), 2) = nodeY(:);
nElements = nelx * nely;
elemNodes = zeros(nElements, 4);
elemDOFs  = zeros(nElements, 8);
elemIndex = zeros(nely, nelx);
elemColRow = zeros(nElements, 2);
num = 0;
for elx = 1:nelx
    for ely = 1:nely
        num = num + 1;
        elemIndex(ely, elx) = num;
        elemColRow(num,:) = [elx, ely];
        n1 = nodeID(ely,     elx);
        n2 = nodeID(ely,     elx + 1);
        n3 = nodeID(ely + 1, elx + 1);
        n4 = nodeID(ely + 1, elx);
        elemNodes(num,:) = [n1, n2, n3, n4];
        elemDOFs(num,:) = [2*n1-1, 2*n1, ...
                           2*n2-1, 2*n2, ...
                           2*n3-1, 2*n3, ...
                           2*n4-1, 2*n4];
    end
end
nDOF = 2 * nNodes;
[XC, YC] = meshgrid(xCenters, yCenters);
rearTri = geometryData.rearTri;
mainQuad = geometryData.mainQuad;
batteryPoly = geometryData.batteryPoly;
anchor = geometryData.anchor;
insideRear = inpolygon(XC, YC, rearTri(:,1), rearTri(:,2));
insideMain = inpolygon(XC, YC, mainQuad(:,1), mainQuad(:,2));
activeMask = insideRear | insideMain;
insideBattery = inpolygon(XC, YC, batteryPoly(:,1), batteryPoly(:,2));
batteryVoidMask = activeMask & insideBattery;
passiveSolidMask = false(size(activeMask));
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, P.RA,     anchor.RA_radius);
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, P.BB,     anchor.BB_radius);
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, P.ST,     anchor.ST_radius);
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, P.HT_top, anchor.HT_radius);
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, P.HT_bot, anchor.HT_radius);
passiveSolidMask = passiveSolidMask | segmentDistanceMask(XC, YC, P.HT_bot, P.HT_top, anchor.HeadTube_halfWidth);
passiveSolidMask = activeMask & passiveSolidMask & ~batteryVoidMask;
outsideDesignMask = ~activeMask;
passiveVoidMask = outsideDesignMask | batteryVoidMask;
maskCode = zeros(size(activeMask));
maskCode(activeMask) = 1;
maskCode(passiveVoidMask) = 0;
maskCode(passiveSolidMask) = 2;
maskVector = maskCode(:);
activeElementIDs = find(maskVector == 1);
passiveVoidElementIDs = find(maskVector == 0);
passiveSolidElementIDs = find(maskVector == 2);
materialElementIDs = find(maskVector ~= 0);
materialNodeIDs = unique(elemNodes(materialElementIDs, :));
batteryCenter = mean(batteryPoly, 1);
supportSpecs = struct([]);
supportSpecs(1).name = 'RA fixed support';
supportSpecs(1).targetPoint = P.RA;
supportSpecs(1).nodeID = findClosestNode(nodeCoords, P.RA, materialNodeIDs);
supportSpecs(1).fixUx = true;
supportSpecs(1).fixUy = true;
supportSpecs(2).name = 'HT_bot roller support';
supportSpecs(2).targetPoint = P.HT_bot;
supportSpecs(2).nodeID = findClosestNode(nodeCoords, P.HT_bot, materialNodeIDs);
supportSpecs(2).fixUx = false;
supportSpecs(2).fixUy = true;
fixedDOFs = [];
for i = 1:numel(supportSpecs)
    nid = supportSpecs(i).nodeID;
    if supportSpecs(i).fixUx
        fixedDOFs(end+1) = 2*nid - 1;
    end
    if supportSpecs(i).fixUy
        fixedDOFs(end+1) = 2*nid;
    end
end
fixedDOFs = unique(fixedDOFs);
allDOFs = 1:nDOF;
freeDOFs = setdiff(allDOFs, fixedDOFs);
loadSpecs = struct([]);
loadSpecs(1).name = 'rider load at ST';
loadSpecs(1).targetPoint = P.ST;
loadSpecs(1).nodeID = findClosestNode(nodeCoords, P.ST, materialNodeIDs);
loadSpecs(1).Fx = loads.rider_Fx_N;
loadSpecs(1).Fy = loads.rider_Fy_N;
loadSpecs(2).name = 'battery weight near down tube';
loadSpecs(2).targetPoint = batteryCenter;
loadSpecs(2).nodeID = findClosestNode(nodeCoords, batteryCenter, materialNodeIDs);
loadSpecs(2).Fx = loads.battery_Fx_N;
loadSpecs(2).Fy = loads.battery_Fy_N;
loadSpecs(3).name = 'motor / pedalling load at BB';
loadSpecs(3).targetPoint = P.BB;
loadSpecs(3).nodeID = findClosestNode(nodeCoords, P.BB, materialNodeIDs);
loadSpecs(3).Fx = loads.motor_Fx_N;
loadSpecs(3).Fy = loads.motor_Fy_N;
F = sparse(nDOF, 1);
for i = 1:numel(loadSpecs)
    nid = loadSpecs(i).nodeID;
    F(2*nid - 1) = F(2*nid - 1) + loadSpecs(i).Fx;
    F(2*nid)     = F(2*nid)     + loadSpecs(i).Fy;
end
fig1 = figure('Color','w','Name','Q4 FE mesh with loads and supports');
hold on; axis equal; box on;
title(sprintf('Structured Q4 FE mesh for e-bike frame: %d x %d elements', nelx, nely));
xlabel('x position (mm)'); ylabel('y position (mm)');
imagesc(xCenters, yCenters, maskCode);
set(gca, 'YDir', 'normal');
colormap(gca, [0.96 0.96 0.96; 0.75 0.85 1.00; 0.20 0.45 0.95]);
cb = colorbar;
cb.Ticks = [0, 1, 2];
cb.TickLabels = {'void/outside','active','passive solid'};
if plotFullGridLines
    nGridLines = numel(xEdges) + numel(yEdges);
    if nGridLines <= maxGridLinesToPlot
        for ix = 1:numel(xEdges)
            plot([xEdges(ix), xEdges(ix)], [yEdges(1), yEdges(end)], '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.20);
        end
        for iy = 1:numel(yEdges)
            plot([xEdges(1), xEdges(end)], [yEdges(iy), yEdges(iy)], '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.20);
        end
    else
        fprintf('Grid line plotting skipped because mesh is very fine (%d grid lines).\n', nGridLines);
    end
end
plotFrameMembers(P);
plot([batteryPoly(:,1); batteryPoly(1,1)], [batteryPoly(:,2); batteryPoly(1,2)], 'r-', 'LineWidth', 2);
for i = 1:numel(supportSpecs)
    nid = supportSpecs(i).nodeID;
    p = nodeCoords(nid,:);
    plot(p(1), p(2), 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    text(p(1) + 8, p(2) - 12, supportSpecs(i).name, 'FontSize', 8, 'Interpreter', 'none');
end
for i = 1:numel(loadSpecs)
    nid = loadSpecs(i).nodeID;
    p = nodeCoords(nid,:);
    plot(p(1), p(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    drawLoadArrow(p, [loadSpecs(i).Fx, loadSpecs(i).Fy]);
    text(p(1) + 8, p(2) + 10, sprintf('%s\nFx=%.0f N, Fy=%.0f N', ...
        loadSpecs(i).name, loadSpecs(i).Fx, loadSpecs(i).Fy), ...
        'FontSize', 8, 'Interpreter', 'none');
end
xlim([xmin, xmax]); ylim([ymin, ymax]);
exportFigure(fig1, fullfile(outDir, '04_fe_mesh_q4_with_loads_supports.png'));
fig2 = figure('Color','w','Name','Q4 element categories');
imagesc(xCenters, yCenters, maskCode);
set(gca, 'YDir', 'normal'); axis equal tight; box on;
colormap(gca, [0.96 0.96 0.96; 0.75 0.85 1.00; 0.20 0.45 0.95]);
cb = colorbar;
cb.Ticks = [0, 1, 2];
cb.TickLabels = {'void/outside','active','passive solid'};
title('Q4 element categories for topology optimization');
xlabel('x position (mm)'); ylabel('y position (mm)');
hold on;
plotFrameMembers(P);
plot([batteryPoly(:,1); batteryPoly(1,1)], [batteryPoly(:,2); batteryPoly(1,2)], 'r-', 'LineWidth', 2);
exportFigure(fig2, fullfile(outDir, '05_q4_element_categories.png'));
if ~isempty(activeElementIDs)
    centerPoint = mean([P.RA; P.BB; P.ST; P.HT_top; P.HT_bot], 1);
    elemCentersVector = [XC(:), YC(:)];
    activeCenters = elemCentersVector(activeElementIDs, :);
    [~, localID] = min(sum((activeCenters - centerPoint).^2, 2));
    exampleElem = activeElementIDs(localID);
else
    exampleElem = 1;
end
fig3 = figure('Color','w','Name','Q4 element numbering example');
hold on; axis equal; box on; grid on;
title(sprintf('Q4 element numbering example: element ID %d', exampleElem));
xlabel('x position (mm)'); ylabel('y position (mm)');
nodesEx = elemNodes(exampleElem,:);
coordsEx = nodeCoords(nodesEx,:);
patch(coordsEx([1 2 3 4],1), coordsEx([1 2 3 4],2), [0.85 0.92 1.0], ...
    'EdgeColor', 'k', 'LineWidth', 2);
plot(coordsEx(:,1), coordsEx(:,2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
localNames = {'n1 bottom-left', 'n2 bottom-right', 'n3 top-right', 'n4 top-left'};
for i = 1:4
    nid = nodesEx(i);
    dofs = [2*nid-1, 2*nid];
    text(coordsEx(i,1), coordsEx(i,2), sprintf('  %s\n  node %d\n  DOF [%d,%d]', ...
        localNames{i}, nid, dofs(1), dofs(2)), ...
        'FontSize', 9, 'Interpreter', 'none');
end
xPad = 2 * (xEdges(2) - xEdges(1));
yPad = 2 * (yEdges(2) - yEdges(1));
xlim([min(coordsEx(:,1))-xPad, max(coordsEx(:,1))+xPad]);
ylim([min(coordsEx(:,2))-yPad, max(coordsEx(:,2))+yPad]);
exportFigure(fig3, fullfile(outDir, '06_q4_element_node_dof_numbering_example.png'));
nodePreview = [(1:min(50,nNodes))', nodeCoords(1:min(50,nNodes),:)];
elemPreview = [(1:min(50,nElements))', elemNodes(1:min(50,nElements),:), elemDOFs(1:min(50,nElements),:), maskVector(1:min(50,nElements))];
writeMatrixCompat(nodePreview, fullfile(outDir, 'node_coordinates_preview_first50.csv'));
writeMatrixCompat(elemPreview, fullfile(outDir, 'element_connectivity_dofs_preview_first50.csv'));
feMesh = struct();
feMesh.units = 'mm, N, MPa';
feMesh.description = 'Structured Q4 mesh generated from MATLAB CAD / geometric model';
feMesh.nelx = nelx;
feMesh.nely = nely;
feMesh.nNodes = nNodes;
feMesh.nElements = nElements;
feMesh.nDOF = nDOF;
feMesh.xEdges = xEdges;
feMesh.yEdges = yEdges;
feMesh.xCenters = xCenters;
feMesh.yCenters = yCenters;
feMesh.XC = XC;
feMesh.YC = YC;
feMesh.nodeID = nodeID;
feMesh.nodeCoords = nodeCoords;
feMesh.elemIndex = elemIndex;
feMesh.elemColRow = elemColRow;
feMesh.elemNodes = elemNodes;
feMesh.elemDOFs = elemDOFs;
feMesh.maskCode = maskCode;
feMesh.maskVector = maskVector;
feMesh.activeElementIDs = activeElementIDs;
feMesh.passiveVoidElementIDs = passiveVoidElementIDs;
feMesh.passiveSolidElementIDs = passiveSolidElementIDs;
feMesh.materialElementIDs = materialElementIDs;
feMesh.materialNodeIDs = materialNodeIDs;
feMesh.exampleElem = exampleElem;
bcData = struct();
bcData.description = 'Starter loads/supports for Step 4 baseline FEM; edit if needed';
bcData.supportSpecs = supportSpecs;
bcData.loadSpecs = loadSpecs;
bcData.F = F;
bcData.fixedDOFs = fixedDOFs;
bcData.freeDOFs = freeDOFs;
bcData.allDOFs = allDOFs;
save(fullfile(outDir, 'ebike_fe_mesh_step3.mat'), 'feMesh', 'bcData', 'geometryData');
fprintf('\nSTEP 3 - Q4 FE mesh generated successfully.\n');
fprintf('Output folder: %s\n', outDir);
fprintf('\nMesh preset: %s\n', meshPreset);
fprintf('nelx x nely      = %d x %d\n', nelx, nely);
fprintf('Total Q4 elements = %d\n', nElements);
fprintf('Total nodes       = %d\n', nNodes);
fprintf('Total DOFs        = %d\n', nDOF);
fprintf('Element size      = %.3f mm x %.3f mm\n', xEdges(2)-xEdges(1), yEdges(2)-yEdges(1));
fprintf('\nElement mask counts:\n');
fprintf('  Passive void / outside / battery cavity = %d\n', numel(passiveVoidElementIDs));
fprintf('  Active design elements                  = %d\n', numel(activeElementIDs));
fprintf('  Passive solid anchor elements           = %d\n', numel(passiveSolidElementIDs));
fprintf('\nExample element 1:\n');
fprintf('  Nodes = [%d %d %d %d]\n', elemNodes(1,1), elemNodes(1,2), elemNodes(1,3), elemNodes(1,4));
fprintf('  DOFs  = [%d %d %d %d %d %d %d %d]\n', elemDOFs(1,:));
fprintf('\nSupport nodes and fixed DOFs:\n');
for i = 1:numel(supportSpecs)
    fprintf('  %-22s node = %d at [%.2f, %.2f], fixUx=%d, fixUy=%d\n', ...
        supportSpecs(i).name, supportSpecs(i).nodeID, ...
        nodeCoords(supportSpecs(i).nodeID,1), nodeCoords(supportSpecs(i).nodeID,2), ...
        supportSpecs(i).fixUx, supportSpecs(i).fixUy);
end
fprintf('  fixedDOFs = '); fprintf('%d ', fixedDOFs); fprintf('\n');
fprintf('\nLoad nodes:\n');
for i = 1:numel(loadSpecs)
    fprintf('  %-30s node = %d at [%.2f, %.2f], Fx=%.1f N, Fy=%.1f N\n', ...
        loadSpecs(i).name, loadSpecs(i).nodeID, ...
        nodeCoords(loadSpecs(i).nodeID,1), nodeCoords(loadSpecs(i).nodeID,2), ...
        loadSpecs(i).Fx, loadSpecs(i).Fy);
end
fprintf('\nFiles exported:\n');
fprintf('  04_fe_mesh_q4_with_loads_supports.png\n');
fprintf('  05_q4_element_categories.png\n');
fprintf('  06_q4_element_node_dof_numbering_example.png\n');
fprintf('  node_coordinates_preview_first50.csv\n');
fprintf('  element_connectivity_dofs_preview_first50.csv\n');
fprintf('  ebike_fe_mesh_step3.mat\n\n');
function mask = pointRadiusMask(X, Y, p, r)
    mask = (X - p(1)).^2 + (Y - p(2)).^2 <= r^2;
end
function mask = segmentDistanceMask(X, Y, p1, p2, halfWidth)
    vx = p2(1) - p1(1);
    vy = p2(2) - p1(2);
    wx = X - p1(1);
    wy = Y - p1(2);
    c1 = wx*vx + wy*vy;
    c2 = vx^2 + vy^2;
    t = max(0, min(1, c1./c2));
    projX = p1(1) + t*vx;
    projY = p1(2) + t*vy;
    dist2 = (X - projX).^2 + (Y - projY).^2;
    mask = dist2 <= halfWidth^2;
end
function node = findClosestNode(nodeCoords, targetPoint, candidateNodeIDs)
    if nargin < 3 || isempty(candidateNodeIDs)
        candidateNodeIDs = (1:size(nodeCoords,1))';
    end
    candidateCoords = nodeCoords(candidateNodeIDs, :);
    d2 = (candidateCoords(:,1) - targetPoint(1)).^2 + (candidateCoords(:,2) - targetPoint(2)).^2;
    [~, idx] = min(d2);
    node = candidateNodeIDs(idx);
end
function plotFrameMembers(P)
    pairs = {P.RA, P.BB; P.RA, P.ST; P.BB, P.ST; P.ST, P.HT_top; P.BB, P.HT_bot; P.HT_bot, P.HT_top};
    for k = 1:size(pairs,1)
        p1 = pairs{k,1};
        p2 = pairs{k,2};
        plot([p1(1), p2(1)], [p1(2), p2(2)], 'k-', 'LineWidth', 1.5);
    end
end
function drawLoadArrow(point, forceVector)
    if norm(forceVector) < eps
        return;
    end
    arrowLength = 55;
    direction = forceVector / norm(forceVector);
    startPt = point - arrowLength * direction;
    quiver(startPt(1), startPt(2), arrowLength*direction(1), arrowLength*direction(2), ...
        0, 'r', 'LineWidth', 2, 'MaxHeadSize', 1.8);
end
function exportFigure(figHandle, filename)
    try
        exportgraphics(figHandle, filename, 'Resolution', 300);
    catch
        saveas(figHandle, filename);
    end
end
function writeMatrixCompat(A, filename)
    try
        writematrix(A, filename);
    catch
        csvwrite(filename, A);
    end
end
