clear; 
clc; 
close all;

volfrac = 0.45;
penal = 3.0;
rmin = 3.0;
move = 0.15;
xmin = 1e-3;
maxIter = 120;
tolx = 0.01;
useBaselineEminRatio = true;
customEminRatio = 1e-3;
thresholdDensity = 0.50;
stressDensityCutoff = 0.30;
showElementEdges = false;
printEvery = 1;
baselineFolder = 'ebike_baseline_fem_outputs';
baselineMatName = 'ebike_baseline_fem_step4_5.mat';
outDir = 'ebike_topology_optimization_outputs';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
baselineFile = fullfile(baselineFolder, baselineMatName);
if ~exist(baselineFile, 'file')
    error(['Cannot find baseline file: %s\n' ...
           'Run baseline_fem.m first, or check that this script is in the same folder as ebike_baseline_fem_outputs.'], baselineFile);
end
load(baselineFile, 'feMesh', 'geometryData', 'baselineData');
nelx = double(feMesh.nelx);
nely = double(feMesh.nely);
nElements = double(feMesh.nElements);
nNodes = double(feMesh.nNodes);
nDOF = double(feMesh.nDOF);
nodeCoords = feMesh.nodeCoords;
elemNodes = double(feMesh.elemNodes);
edofMat = double(feMesh.elemDOFs);
xCenters = feMesh.xCenters;
yCenters = feMesh.yCenters;
xEdges = feMesh.xEdges;
yEdges = feMesh.yEdges;
dx = xEdges(2) - xEdges(1);
dy = yEdges(2) - yEdges(1);
elementArea = dx * dy;
P = geometryData.points;
batteryPoly = geometryData.batteryPoly;
F = baselineData.F;
if ~issparse(F)
    F = sparse(F);
end
F = F(:);
fixedDOFs = double(baselineData.fixedDOFs(:));
freeDOFs = double(baselineData.freeDOFs(:));
KE0 = baselineData.KE0;
D = baselineData.D;
material = baselineData.material;
thickness = double(material.thickness_mm);
if useBaselineEminRatio && isfield(baselineData, 'EminRatio')
    EminRatio = double(baselineData.EminRatio);
else
    EminRatio = customEminRatio;
end
if isfield(baselineData, 'maskVectorClean')
    maskVector = double(baselineData.maskVectorClean(:));
    maskCode = double(baselineData.maskCodeClean);
else
    maskVector = double(feMesh.maskVector(:));
    maskCode = double(feMesh.maskCode);
end
activeElementIDs = find(maskVector == 1);
passiveSolidElementIDs = find(maskVector == 2);
passiveVoidElementIDs = find(maskVector == 0);
activeMaskGrid = reshape(maskVector == 1, nely, nelx);
passiveSolidMaskGrid = reshape(maskVector == 2, nely, nelx);
passiveVoidMaskGrid = reshape(maskVector == 0, nely, nelx);
nActive = numel(activeElementIDs);
nPassiveSolid = numel(passiveSolidElementIDs);
nPassiveVoid = numel(passiveVoidElementIDs);
if nActive == 0
    error('No active design elements were found. Check Step 3/4 masks.');
end
if volfrac <= xmin || volfrac > 1
    error('volfrac must be between xmin and 1. Current volfrac = %.4f, xmin = %.4f.', volfrac, xmin);
end
baselineCompliance = full(baselineData.compliance);
baselineMaxDisp = full(baselineData.maxDispMaterial);
baselineMaterialVolume = full(baselineData.materialVolumeBaseline);
iK = reshape(kron(edofMat, ones(8,1))', 64*nElements, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nElements, 1);
[H, Hs] = buildActiveDensityFilter(activeElementIDs, maskVector, nelx, nely, rmin);
x = volfrac * ones(nActive, 1);
xPhys = makeFullDensityVector(x, H, Hs, nElements, activeElementIDs, passiveSolidElementIDs);
history.iter = zeros(maxIter, 1);
history.compliance = zeros(maxIter, 1);
history.strainEnergy = zeros(maxIter, 1);
history.meanActiveDensity = zeros(maxIter, 1);
history.activeVolumeFraction = zeros(maxIter, 1);
history.totalMaterialVolume_mm3 = zeros(maxIter, 1);
history.maxMaterialDisp_mm = zeros(maxIter, 1);
history.change = zeros(maxIter, 1);
history.physicalChange = zeros(maxIter, 1);
fprintf('\nSTEP 6 - TOPOLOGY OPTIMIZATION STARTED\n');
fprintf('Input baseline file: %s\n', baselineFile);
fprintf('Output folder:      %s\n\n', outDir);
fprintf('Optimization settings:\n');
fprintf('  volfrac active domain = %.3f\n', volfrac);
fprintf('  penal                 = %.3f\n', penal);
fprintf('  rmin                  = %.3f elements\n', rmin);
fprintf('  move                  = %.3f\n', move);
fprintf('  xmin                  = %.4g\n', xmin);
fprintf('  maxIter               = %d\n', maxIter);
fprintf('  tolx                  = %.4f\n', tolx);
fprintf('  EminRatio             = %.4g\n\n', EminRatio);
fprintf('Design domain summary:\n');
fprintf('  nelx x nely           = %d x %d\n', nelx, nely);
fprintf('  total elements        = %d\n', nElements);
fprintf('  active design elems   = %d\n', nActive);
fprintf('  passive solid elems   = %d\n', nPassiveSolid);
fprintf('  passive void elems    = %d\n\n', nPassiveVoid);
fprintf('Iter | Compliance (Nmm) | Active Vol | Max Disp (mm) | Change\n');
fprintf('\n');
U = zeros(nDOF, 1);
lastXPhysActive = xPhys(activeElementIDs);
for loop = 1:maxIter
    xPhys = makeFullDensityVector(x, H, Hs, nElements, activeElementIDs, passiveSolidElementIDs);
    xPhysActive = xPhys(activeElementIDs);
    [U, compliance, elementComplianceFullMaterial, stiffnessScale] = solveFEMTopology( ...
        xPhys, penal, EminRatio, KE0, edofMat, iK, jK, nDOF, freeDOFs, F);
    [maxDispMaterial, ~, ~] = computeMaxMaterialDisplacement(U, elemNodes, xPhys, passiveSolidElementIDs, stressDensityCutoff);
    dcPhys = zeros(nElements, 1);
    dcPhys(activeElementIDs) = -penal * (1 - EminRatio) .* ...
        (max(xPhys(activeElementIDs), xmin).^(penal - 1)) .* elementComplianceFullMaterial(activeElementIDs);
    dvPhysActive = ones(nActive, 1);
    dc = H' * (dcPhys(activeElementIDs) ./ Hs);
    dv = H' * (dvPhysActive ./ Hs);
    xOld = x;
    x = OCUpdateDensityFilter(xOld, dc, dv, H, Hs, volfrac, move, xmin);
    xPhysNew = makeFullDensityVector(x, H, Hs, nElements, activeElementIDs, passiveSolidElementIDs);
    xPhysActiveNew = xPhysNew(activeElementIDs);
    change = max(abs(x(:) - xOld(:)));
    physicalChange = max(abs(xPhysActiveNew(:) - lastXPhysActive(:)));
    lastXPhysActive = xPhysActiveNew;
    activeVolFrac = mean(xPhysActiveNew);
    totalMaterialVolume = (sum(xPhysActiveNew) + nPassiveSolid) * elementArea * thickness;
    history.iter(loop) = loop;
    history.compliance(loop) = full(compliance);
    history.strainEnergy(loop) = 0.5 * full(compliance);
    history.meanActiveDensity(loop) = mean(xPhysActiveNew);
    history.activeVolumeFraction(loop) = activeVolFrac;
    history.totalMaterialVolume_mm3(loop) = totalMaterialVolume;
    history.maxMaterialDisp_mm(loop) = maxDispMaterial;
    history.change(loop) = change;
    history.physicalChange(loop) = physicalChange;
    if mod(loop, printEvery) == 0 || loop == 1
        fprintf('%4d | %16.6e | %10.4f | %13.6e | %8.4f\n', ...
            loop, full(compliance), activeVolFrac, maxDispMaterial, change);
    end
    if change < tolx
        fprintf('Converged because max design change %.4f < tolx %.4f.\n', change, tolx);
        break;
    end
end
nIter = loop;
fields = fieldnames(history);
for k = 1:numel(fields)
    history.(fields{k}) = history.(fields{k})(1:nIter);
end
xPhysFinal = makeFullDensityVector(x, H, Hs, nElements, activeElementIDs, passiveSolidElementIDs);
[UFinal, finalCompliance, elementComplianceFullMaterialFinal, stiffnessScaleFinal, KFinal] = solveFEMTopology( ...
    xPhysFinal, penal, EminRatio, KE0, edofMat, iK, jK, nDOF, freeDOFs, F);
finalStrainEnergy = 0.5 * full(finalCompliance);
Ux = UFinal(1:2:end);
Uy = UFinal(2:2:end);
dispMag = sqrt(Ux.^2 + Uy.^2);
[maxDispMaterial, maxDispNode, materialNodeIDsFinal] = computeMaxMaterialDisplacement( ...
    UFinal, elemNodes, xPhysFinal, passiveSolidElementIDs, stressDensityCutoff);
maxDispGlobal = max(dispMag);
UeFinal = UFinal(edofMat);
elementComplianceActualFinal = stiffnessScaleFinal .* elementComplianceFullMaterialFinal;
elementStrainEnergyFinal = 0.5 * elementComplianceActualFinal;
coordsRect = [0, 0; dx, 0; dx, dy; 0, dy];
[~, Bcenter] = q4CenterBMatrix(coordsRect);
stressCenterFinal = zeros(nElements, 3);
vonMisesFinal = zeros(nElements, 1);
stressPlotElementIDs = find((xPhysFinal >= stressDensityCutoff) | (maskVector == 2));
for kk = 1:numel(stressPlotElementIDs)
    e = stressPlotElementIDs(kk);
    ue = UeFinal(e, :)';
    sigma = D * Bcenter * ue;
    stressCenterFinal(e, :) = sigma';
    sx = sigma(1); sy = sigma(2); txy = sigma(3);
    vonMisesFinal(e) = sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2);
end
if isempty(stressPlotElementIDs)
    maxVonMises = NaN;
else
    maxVonMises = max(vonMisesFinal(stressPlotElementIDs));
end
activePhysicalVolume = sum(xPhysFinal(activeElementIDs)) * elementArea * thickness;
passiveSolidVolume = nPassiveSolid * elementArea * thickness;
finalMaterialVolume = activePhysicalVolume + passiveSolidVolume;
finalVolumeRatioToBaseline = finalMaterialVolume / baselineMaterialVolume;
finalVolumeReductionPercent = 100 * (1 - finalVolumeRatioToBaseline);
binarySolidActive = xPhysFinal(activeElementIDs) >= thresholdDensity;
binaryMaterialVolume = (nnz(binarySolidActive) + nPassiveSolid) * elementArea * thickness;
edgeColor = 'none';
if showElementEdges
    edgeColor = [0.75 0.75 0.75];
end
xPhysGrid = reshape(xPhysFinal, nely, nelx);
maskMaterial = reshape(maskVector ~= 0, nely, nelx);
binaryGrid = zeros(nely, nelx);
binaryGrid(activeMaskGrid) = xPhysGrid(activeMaskGrid) >= thresholdDensity;
binaryGrid(passiveSolidMaskGrid) = 1;
fig12 = figure('Color','w','Name','Step 6 optimization formulation');
hold on; axis equal; box on;
imagesc(xCenters, yCenters, maskCode);
set(gca, 'YDir', 'normal');
colormap(gca, [0.96 0.96 0.96; 0.72 0.86 1.00; 0.20 0.45 0.95]);
cb = colorbar;
cb.Ticks = [0, 1, 2];
cb.TickLabels = {'passive void','active design','passive solid'};
plotFrameMembers(P);
plot([batteryPoly(:,1); batteryPoly(1,1)], [batteryPoly(:,2); batteryPoly(1,2)], 'r-', 'LineWidth', 2);
plotLoadsSupports(nodeCoords, baselineData.supportSpecs, baselineData.loadSpecs);
title('Step 6 formulation: active domain, passive regions, loads and supports');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig12, fullfile(outDir, '12_topopt_formulation_mask.png'));
fig13 = figure('Color','w','Name','Topology optimization convergence history');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(history.iter, history.compliance, '-o', 'LineWidth', 1.2, 'MarkerSize', 3);
grid on; xlabel('iteration'); ylabel('C = F^T U (N mm)');
title('Compliance history');
nexttile;
plot(history.iter, history.activeVolumeFraction, '-o', 'LineWidth', 1.2, 'MarkerSize', 3);
hold on; yline(volfrac, '--', 'target volfrac');
grid on; xlabel('iteration'); ylabel('active volume fraction');
title('Volume constraint history');
nexttile;
plot(history.iter, history.change, '-o', 'LineWidth', 1.2, 'MarkerSize', 3);
hold on; yline(tolx, '--', 'tolx');
grid on; xlabel('iteration'); ylabel('max change');
title('Design variable change history');
exportFigure(fig13, fullfile(outDir, '13_topopt_iteration_history.png'));
fig14 = figure('Color','w','Name','Optimized density map');
hold on; axis equal; box on;
imagesc(xCenters, yCenters, xPhysGrid, [0 1]);
set(gca, 'YDir', 'normal');
colormap(gca, gray);
colorbar;
contour(xCenters, yCenters, double(maskMaterial), [0.5 0.5], 'k-', 'LineWidth', 1.0);
plotFrameMembers(P);
plot([batteryPoly(:,1); batteryPoly(1,1)], [batteryPoly(:,2); batteryPoly(1,2)], 'r-', 'LineWidth', 2);
plotLoadsSupports(nodeCoords, baselineData.supportSpecs, baselineData.loadSpecs);
title(sprintf('Optimized density map, volfrac = %.2f, penal = %.1f', volfrac, penal));
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig14, fullfile(outDir, '14_optimized_density_map.png'));
fig15 = figure('Color','w','Name','Optimized binary topology threshold');
hold on; axis equal; box on;
imagesc(xCenters, yCenters, binaryGrid, [0 1]);
set(gca, 'YDir', 'normal');
colormap(gca, gray);
colorbar;
contour(xCenters, yCenters, double(maskMaterial), [0.5 0.5], 'k-', 'LineWidth', 1.0);
plotFrameMembers(P);
plot([batteryPoly(:,1); batteryPoly(1,1)], [batteryPoly(:,2); batteryPoly(1,2)], 'r-', 'LineWidth', 2);
plotLoadsSupports(nodeCoords, baselineData.supportSpecs, baselineData.loadSpecs);
title(sprintf('Binary redesign guide, density threshold = %.2f', thresholdDensity));
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig15, fullfile(outDir, '15_optimized_topology_threshold.png'));
plotElementIDs = find((xPhysFinal >= stressDensityCutoff) | (maskVector == 2));
if isempty(plotElementIDs)
    plotElementIDs = passiveSolidElementIDs;
end
if maxDispMaterial > eps
    modelSize = max([max(xEdges)-min(xEdges), max(yEdges)-min(yEdges)]);
    defScale = 0.10 * modelSize / maxDispMaterial;
else
    defScale = 1.0;
end
deformedCoords = nodeCoords + defScale * [Ux, Uy];
fig16 = figure('Color','w','Name','Optimized deformed shape');
hold on; axis equal; box on;
patch('Faces', elemNodes(plotElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceColor', [0.90 0.90 0.90], 'EdgeColor', [0.80 0.80 0.80], 'LineWidth', 0.15);
patch('Faces', elemNodes(plotElementIDs,:), 'Vertices', deformedCoords, ...
      'FaceVertexCData', dispMag, 'FaceColor', 'interp', 'EdgeColor', edgeColor);
colorbar;
plotFrameMembers(P);
title(sprintf('Optimized deformed shape, scale = %.1f x', defScale));
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig16, fullfile(outDir, '16_optimized_deformed_shape.png'));
fig17 = figure('Color','w','Name','Optimized displacement magnitude');
patch('Faces', elemNodes(plotElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', dispMag, 'FaceColor', 'interp', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotLoadsSupports(nodeCoords, baselineData.supportSpecs, baselineData.loadSpecs);
title('Optimized topology displacement magnitude');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig17, fullfile(outDir, '17_optimized_displacement_magnitude.png'));
fig18 = figure('Color','w','Name','Optimized strain energy');
patch('Faces', elemNodes(plotElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', elementStrainEnergyFinal(plotElementIDs), ...
      'FaceColor', 'flat', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotFrameMembers(P);
title('Optimized topology element strain energy');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig18, fullfile(outDir, '18_optimized_element_strain_energy.png'));
fig19 = figure('Color','w','Name','Optimized von Mises stress');
patch('Faces', elemNodes(stressPlotElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', vonMisesFinal(stressPlotElementIDs), ...
      'FaceColor', 'flat', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotFrameMembers(P);
title(sprintf('Optimized topology von Mises stress, density cutoff = %.2f', stressDensityCutoff));
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig19, fullfile(outDir, '19_optimized_von_mises_center.png'));
fig20 = figure('Color','w','Name','Baseline vs optimized comparison');
metrics = [finalCompliance / baselineCompliance, ...
           maxDispMaterial / baselineMaxDisp, ...
           finalMaterialVolume / baselineMaterialVolume];
bar(metrics);
set(gca, 'XTickLabel', {'Compliance ratio','Max disp. ratio','Material volume ratio'});
yline(1, '--', 'baseline');
grid on; ylabel('optimized / baseline');
title('Baseline solid domain vs optimized topology');
exportFigure(fig20, fullfile(outDir, '20_baseline_vs_optimized_comparison.png'));
historyTable = table(history.iter, history.compliance, history.strainEnergy, ...
    history.meanActiveDensity, history.activeVolumeFraction, history.totalMaterialVolume_mm3, ...
    history.maxMaterialDisp_mm, history.change, history.physicalChange, ...
    'VariableNames', {'iteration','compliance_Nmm','strain_energy_Nmm','mean_active_density', ...
    'active_volume_fraction','total_material_volume_mm3','max_material_disp_mm','design_change','physical_density_change'});
writetable(historyTable, fullfile(outDir, 'topopt_iteration_history.csv'));
summaryNames = {'volfrac_active_domain'; 'penal'; 'rmin_elements'; 'move_limit'; 'xmin'; ...
                'iterations'; 'baseline_compliance_Nmm'; 'optimized_compliance_Nmm'; ...
                'compliance_ratio_opt_over_baseline'; 'baseline_max_disp_mm'; 'optimized_max_disp_mm'; ...
                'disp_ratio_opt_over_baseline'; 'baseline_material_volume_mm3'; ...
                'optimized_physical_material_volume_mm3'; 'material_volume_ratio_opt_over_baseline'; ...
                'material_volume_reduction_percent'; 'binary_material_volume_mm3_thresholded'; ...
                'max_von_mises_MPa_qualitative'};
summaryValues = [volfrac; penal; rmin; move; xmin; nIter; baselineCompliance; full(finalCompliance); ...
                 full(finalCompliance)/baselineCompliance; baselineMaxDisp; maxDispMaterial; ...
                 maxDispMaterial/baselineMaxDisp; baselineMaterialVolume; finalMaterialVolume; ...
                 finalVolumeRatioToBaseline; finalVolumeReductionPercent; binaryMaterialVolume; maxVonMises];
summaryTable = table(summaryNames, summaryValues, 'VariableNames', {'item','value'});
writetable(summaryTable, fullfile(outDir, 'topopt_summary.csv'));
parameterNames = {'material_case'; 'E_MPa'; 'nu'; 'thickness_mm'; 'EminRatio'; 'filter_type'; 'optimizer'};
parameterValues = {char(material.name); ...
                   sprintf('%.6g', material.E_MPa); ...
                   sprintf('%.6g', material.nu); ...
                   sprintf('%.6g', material.thickness_mm); ...
                   sprintf('%.6g', EminRatio); ...
                   'active-domain density filter'; ...
                   'OC update'};
parameterTable = table(parameterNames, parameterValues, 'VariableNames', {'parameter','value'});
writetable(parameterTable, fullfile(outDir, 'topopt_parameters.csv'));
topOptData = struct();
topOptData.settings = struct('volfrac', volfrac, 'penal', penal, 'rmin', rmin, 'move', move, ...
    'xmin', xmin, 'maxIter', maxIter, 'tolx', tolx, 'EminRatio', EminRatio, ...
    'thresholdDensity', thresholdDensity, 'stressDensityCutoff', stressDensityCutoff);
topOptData.x = x;
topOptData.xPhys = xPhysFinal;
topOptData.xPhysGrid = xPhysGrid;
topOptData.binaryGrid = binaryGrid;
topOptData.U = UFinal;
topOptData.compliance = full(finalCompliance);
topOptData.strainEnergyTotal = finalStrainEnergy;
topOptData.maxDispMaterial = maxDispMaterial;
topOptData.maxDispNode = maxDispNode;
topOptData.maxDispGlobal = maxDispGlobal;
topOptData.elementComplianceFullMaterial = elementComplianceFullMaterialFinal;
topOptData.elementStrainEnergy = elementStrainEnergyFinal;
topOptData.stiffnessScale = stiffnessScaleFinal;
topOptData.stressCenter = stressCenterFinal;
topOptData.vonMises = vonMisesFinal;
topOptData.maxVonMises = maxVonMises;
topOptData.finalMaterialVolume = finalMaterialVolume;
topOptData.finalVolumeRatioToBaseline = finalVolumeRatioToBaseline;
topOptData.finalVolumeReductionPercent = finalVolumeReductionPercent;
topOptData.binaryMaterialVolume = binaryMaterialVolume;
topOptData.history = history;
topOptData.activeElementIDs = activeElementIDs;
topOptData.passiveSolidElementIDs = passiveSolidElementIDs;
topOptData.passiveVoidElementIDs = passiveVoidElementIDs;
topOptData.materialNodeIDsFinal = materialNodeIDsFinal;
topOptData.filter.H = H;
topOptData.filter.Hs = Hs;
save(fullfile(outDir, 'ebike_topology_optimization_step6.mat'), ...
    'topOptData', 'feMesh', 'geometryData', 'baselineData', '-v7.3');
fprintf('\nSTEP 6 completed successfully.\n');
fprintf('Output folder: %s\n\n', outDir);
fprintf('Final topology optimization summary:\n');
fprintf('  Iterations                         = %d\n', nIter);
fprintf('  Baseline compliance                = %.6e N mm\n', baselineCompliance);
fprintf('  Optimized compliance               = %.6e N mm\n', full(finalCompliance));
fprintf('  Compliance ratio opt/baseline      = %.4f\n', full(finalCompliance)/baselineCompliance);
fprintf('  Baseline max material displacement = %.6e mm\n', baselineMaxDisp);
fprintf('  Optimized max material displacement= %.6e mm\n', maxDispMaterial);
fprintf('  Displacement ratio opt/baseline    = %.4f\n', maxDispMaterial/baselineMaxDisp);
fprintf('  Baseline material volume           = %.6e mm^3\n', baselineMaterialVolume);
fprintf('  Optimized physical material volume = %.6e mm^3\n', finalMaterialVolume);
fprintf('  Material volume ratio opt/baseline = %.4f\n', finalVolumeRatioToBaseline);
fprintf('  Material volume reduction          = %.2f %%\n', finalVolumeReductionPercent);
fprintf('  Qualitative max von Mises          = %.6e MPa\n', maxVonMises);
fprintf('\nFiles exported:\n');
fprintf('  12_topopt_formulation_mask.png\n');
fprintf('  13_topopt_iteration_history.png\n');
fprintf('  14_optimized_density_map.png\n');
fprintf('  15_optimized_topology_threshold.png\n');
fprintf('  16_optimized_deformed_shape.png\n');
fprintf('  17_optimized_displacement_magnitude.png\n');
fprintf('  18_optimized_element_strain_energy.png\n');
fprintf('  19_optimized_von_mises_center.png\n');
fprintf('  20_baseline_vs_optimized_comparison.png\n');
fprintf('  topopt_iteration_history.csv\n');
fprintf('  topopt_summary.csv\n');
fprintf('  topopt_parameters.csv\n');
fprintf('  ebike_topology_optimization_step6.mat\n');
function [H, Hs] = buildActiveDensityFilter(activeElementIDs, maskVector, nelx, nely, rmin)
    nElements = numel(maskVector);
    nActive = numel(activeElementIDs);
    activeMap = zeros(nElements, 1);
    activeMap(activeElementIDs) = 1:nActive;
    maxNeighbors = nActive * (2*ceil(rmin)-1)^2;
    iH = zeros(maxNeighbors, 1);
    jH = zeros(maxNeighbors, 1);
    sH = zeros(maxNeighbors, 1);
    k = 0;
    rceil = ceil(rmin) - 1;
    for ii = 1:nActive
        e = activeElementIDs(ii);
        [row_i, col_i] = ind2sub([nely, nelx], e);
        colMin = max(col_i - rceil, 1);
        colMax = min(col_i + rceil, nelx);
        rowMin = max(row_i - rceil, 1);
        rowMax = min(row_i + rceil, nely);
        for col_j = colMin:colMax
            for row_j = rowMin:rowMax
                ej = sub2ind([nely, nelx], row_j, col_j);
                jj = activeMap(ej);
                if jj > 0
                    fac = rmin - sqrt((double(col_i)-double(col_j))^2 + (double(row_i)-double(row_j))^2);
                    if fac > 0
                        k = k + 1;
                        iH(k) = ii;
                        jH(k) = jj;
                        sH(k) = fac;
                    end
                end
            end
        end
    end
    iH = iH(1:k);
    jH = jH(1:k);
    sH = sH(1:k);
    H = sparse(iH, jH, sH, nActive, nActive);
    Hs = full(sum(H, 2));
    Hs(Hs == 0) = 1;
end
function xPhys = makeFullDensityVector(xActive, H, Hs, nElements, activeElementIDs, passiveSolidElementIDs)
    xPhys = zeros(nElements, 1);
    xPhys(activeElementIDs) = full(H * xActive ./ Hs);
    xPhys(passiveSolidElementIDs) = 1.0;
end
function xNew = OCUpdateDensityFilter(xOld, dc, dv, H, Hs, volfrac, move, xmin)
    l1 = 0;
    l2 = 1e9;
    xNew = xOld;
    for iter = 1:100
        lmid = 0.5 * (l1 + l2);
        updateFactor = sqrt(max(0, -dc ./ max(1e-30, dv * lmid)));
        candidate = xOld .* updateFactor;
        candidate = min(xOld + move, candidate);
        candidate = max(xOld - move, candidate);
        candidate = min(1.0, candidate);
        candidate = max(xmin, candidate);
        candidatePhys = full(H * candidate ./ Hs);
        if mean(candidatePhys) > volfrac
            l1 = lmid;
        else
            l2 = lmid;
        end
        xNew = candidate;
        if (l2 - l1) / max(1, (l1 + l2)) < 1e-4
            break;
        end
    end
end
function [U, compliance, elementComplianceFullMaterial, stiffnessScale, K] = solveFEMTopology( ...
    xPhys, penal, EminRatio, KE0, edofMat, iK, jK, nDOF, freeDOFs, F)
    nElements = numel(xPhys);
    stiffnessScale = EminRatio + (xPhys(:).^penal) * (1 - EminRatio);
    sK = reshape(KE0(:) * stiffnessScale(:)', 64*nElements, 1);
    K = sparse(iK, jK, sK, nDOF, nDOF);
    K = (K + K') / 2;
    U = zeros(nDOF, 1);
    U(freeDOFs) = K(freeDOFs, freeDOFs) \ F(freeDOFs);
    if any(~isfinite(U))
        error('FEM solution contains NaN or Inf. Try increasing EminRatio or checking supports.');
    end
    compliance = full(F' * U);
    Ue = U(edofMat);
    elementComplianceFullMaterial = sum((Ue * KE0) .* Ue, 2);
end
function [maxDispMaterial, maxDispNode, materialNodeIDs] = computeMaxMaterialDisplacement(U, elemNodes, xPhys, passiveSolidElementIDs, densityCutoff)
    Ux = U(1:2:end);
    Uy = U(2:2:end);
    dispMag = sqrt(Ux.^2 + Uy.^2);
    materialElementIDs = find(xPhys >= densityCutoff);
    materialElementIDs = unique([materialElementIDs(:); passiveSolidElementIDs(:)]);
    if isempty(materialElementIDs)
        materialElementIDs = passiveSolidElementIDs(:);
    end
    materialNodeIDs = unique(elemNodes(materialElementIDs, :));
    [maxDispMaterial, localID] = max(dispMag(materialNodeIDs));
    maxDispNode = materialNodeIDs(localID);
end
function [N, B] = q4CenterBMatrix(coords)
    xi = 0;
    eta = 0;
    N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dN_dxi = 0.25 * [-(1-eta),  (1-eta),  (1+eta), -(1+eta)];
    dN_deta = 0.25 * [-(1-xi), -(1+xi),  (1+xi),  (1-xi)];
    J = [dN_dxi; dN_deta] * coords;
    dN_global = J \ [dN_dxi; dN_deta];
    dN_dx = dN_global(1,:);
    dN_dy = dN_global(2,:);
    B = zeros(3, 8);
    for a = 1:4
        B(1, 2*a-1) = dN_dx(a);
        B(2, 2*a)   = dN_dy(a);
        B(3, 2*a-1) = dN_dy(a);
        B(3, 2*a)   = dN_dx(a);
    end
end
function plotFrameMembers(P)
    members = {P.RA, P.BB; P.RA, P.ST; P.BB, P.ST; P.ST, P.HT_top; P.BB, P.HT_bot; P.HT_bot, P.HT_top};
    for i = 1:size(members, 1)
        a = members{i,1}; b = members{i,2};
        plot([a(1), b(1)], [a(2), b(2)], 'k--', 'LineWidth', 1.0);
    end
    keyNames = {'RA','BB','ST','HT top','HT bot'};
    keyPts = [P.RA; P.BB; P.ST; P.HT_top; P.HT_bot];
    plot(keyPts(:,1), keyPts(:,2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
    for i = 1:size(keyPts,1)
        text(keyPts(i,1)+8, keyPts(i,2)+8, keyNames{i}, 'FontSize', 8, 'Interpreter', 'none');
    end
end
function plotLoadsSupports(nodeCoords, supportSpecs, loadSpecs)
    for i = 1:numel(supportSpecs)
        nid = supportSpecs(i).nodeID;
        xy = nodeCoords(nid, :);
        plot(xy(1), xy(2), 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
        shortName = supportSpecs(i).name;
        shortName = strrep(shortName, ' support', '');
        text(xy(1)+12, xy(2)-12, shortName, 'FontSize', 8, 'Interpreter', 'none');
    end
    allFx = [];
    allFy = [];
    for i = 1:numel(loadSpecs)
        allFx(end+1) = loadSpecs(i).Fx;
        allFy(end+1) = loadSpecs(i).Fy;
    end
    maxLoad = max(1, max(sqrt(allFx.^2 + allFy.^2)));
    arrowScale = 90 / maxLoad;
    for i = 1:numel(loadSpecs)
        nid = loadSpecs(i).nodeID;
        xy = nodeCoords(nid, :);
        qx = loadSpecs(i).Fx * arrowScale;
        qy = loadSpecs(i).Fy * arrowScale;
        quiver(xy(1), xy(2), qx, qy, 0, 'r', 'LineWidth', 1.8, 'MaxHeadSize', 0.8);
        plot(xy(1), xy(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
        txt = loadSpecs(i).name;
        txt = strrep(txt, ' load', '');
        text(xy(1)+12, xy(2)+12, txt, 'FontSize', 8, 'Color', 'r', 'Interpreter', 'none');
    end
end
function exportFigure(figHandle, filename)
    set(figHandle, 'PaperPositionMode', 'auto');
    try
        exportgraphics(figHandle, filename, 'Resolution', 300);
    catch
        saveas(figHandle, filename);
    end
end
