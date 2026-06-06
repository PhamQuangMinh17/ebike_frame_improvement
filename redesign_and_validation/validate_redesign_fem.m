clear; clc; close all;
inputFile = fullfile('ebike_redesign_outputs', 'ebike_redesigned_geometry_step7.mat');
outDir = 'ebike_validation_outputs';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
useBaselineEminRatio = true;
customEminRatio = 1e-6;
useBaselinePenal = true;
customPenal = 3.0;
solidThreshold = 0.50;
reprojectLoadsAndSupportsToRedesign = true;
showElementEdges = false;
deformedShapeAutoScale = true;
manualDeformationScale = 30;
exportResolution = 300;
if ~exist(inputFile, 'file')
    error(['Cannot find final clean redesign file:\n  %s\n\n', ...
           'Run interpret_topology_redesign.m first, or check the folder name.'], inputFile);
end
S = load(inputFile);
requiredVars = {'redesignData','feMesh','geometryData','baselineData','topOptData'};
for i = 1:numel(requiredVars)
    if ~isfield(S, requiredVars{i})
        error('Input MAT file does not contain required variable: %s', requiredVars{i});
    end
end
redesignData = S.redesignData;
feMesh = S.feMesh;
geometryData = S.geometryData;
baselineData = S.baselineData;
topOptData = S.topOptData;
nodeCoords = feMesh.nodeCoords;
elemNodes  = double(feMesh.elemNodes);
edofMat    = double(feMesh.elemDOFs);
nElements  = double(feMesh.nElements);
nNodes     = double(feMesh.nNodes);
nDOF       = double(feMesh.nDOF);
nelx       = double(feMesh.nelx);
nely       = double(feMesh.nely);
xCenters   = feMesh.xCenters(:);
yCenters   = feMesh.yCenters(:);
xEdges     = feMesh.xEdges(:);
yEdges     = feMesh.yEdges(:);
dx = xEdges(2) - xEdges(1);
dy = yEdges(2) - yEdges(1);
redesignMaskGrid = orientGrid(logical(redesignData.redesignMask), nelx, nely, 'redesignData.redesignMask');
maskCode = orientGrid(double(baselineData.maskCodeClean), nelx, nely, 'baselineData.maskCodeClean');
xPhysRedesign = double(redesignData.xPhysRedesignVector(:));
if numel(xPhysRedesign) ~= nElements
    error('redesignData.xPhysRedesignVector has %d entries, expected %d.', numel(xPhysRedesign), nElements);
end
solidElementIDs = find(xPhysRedesign >= solidThreshold);
solidNodeIDs = unique(elemNodes(solidElementIDs, :));
if isempty(solidElementIDs)
    error('No solid redesigned elements found. Check the final clean redesign output.');
end
material = baselineData.material;
D = baselineData.D;
KE0 = baselineData.KE0;
if useBaselineEminRatio && isfield(baselineData, 'EminRatio')
    EminRatio = double(baselineData.EminRatio);
else
    EminRatio = customEminRatio;
end
if useBaselinePenal && isfield(baselineData, 'penal')
    penal = double(baselineData.penal);
else
    penal = customPenal;
end
elementArea = double(baselineData.elementArea);
elementVolume = double(baselineData.elementVolume);
if isfield(redesignData, 'redesignVolume_mm3')
    redesignedVolume = double(redesignData.redesignVolume_mm3);
else
    redesignedVolume = sum(xPhysRedesign >= solidThreshold) * elementVolume;
end
baselineVolume = double(baselineData.materialVolumeBaseline);
topologyVolume = double(topOptData.finalMaterialVolume);
P = geometryData.points;
batteryPoly = ensureNx2(geometryData.batteryPoly, 'geometryData.batteryPoly');
if reprojectLoadsAndSupportsToRedesign
    candidateNodesForBC = solidNodeIDs(:);
else
    candidateNodesForBC = baselineData.materialNodeIDs(:);
end
supportSpecsRedesign = baselineData.supportSpecs;
for i = 1:numel(supportSpecsRedesign)
    target = rowPoint(supportSpecsRedesign(i).targetPoint);
    oldNode = supportSpecsRedesign(i).nodeID;
    supportSpecsRedesign(i).originalNodeID = oldNode;
    supportSpecsRedesign(i).nodeID = findClosestNode(nodeCoords, target, candidateNodesForBC);
    supportSpecsRedesign(i).snapDistance_mm = norm(nodeCoords(supportSpecsRedesign(i).nodeID,:) - target);
end
fixedDOFs = [];
for i = 1:numel(supportSpecsRedesign)
    nid = supportSpecsRedesign(i).nodeID;
    if supportSpecsRedesign(i).fixUx
        fixedDOFs(end+1) = 2*nid - 1;
    end
    if supportSpecsRedesign(i).fixUy
        fixedDOFs(end+1) = 2*nid;
    end
end
fixedDOFs = unique(fixedDOFs(:));
freeDOFs = setdiff((1:nDOF)', fixedDOFs);
loadSpecsRedesign = baselineData.loadSpecs;
F = sparse(nDOF, 1);
for i = 1:numel(loadSpecsRedesign)
    target = rowPoint(loadSpecsRedesign(i).targetPoint);
    oldNode = loadSpecsRedesign(i).nodeID;
    loadSpecsRedesign(i).originalNodeID = oldNode;
    loadSpecsRedesign(i).nodeID = findClosestNode(nodeCoords, target, candidateNodesForBC);
    loadSpecsRedesign(i).snapDistance_mm = norm(nodeCoords(loadSpecsRedesign(i).nodeID,:) - target);
    nid = loadSpecsRedesign(i).nodeID;
    F(2*nid - 1) = F(2*nid - 1) + loadSpecsRedesign(i).Fx;
    F(2*nid)     = F(2*nid)     + loadSpecsRedesign(i).Fy;
end
edgeColor = 'none';
if showElementEdges
    edgeColor = [0.70 0.70 0.70];
end
resultRedesign = solveStaticQ4Case(xPhysRedesign, solidElementIDs, solidNodeIDs, ...
    KE0, D, edofMat, elemNodes, nodeCoords, nDOF, fixedDOFs, freeDOFs, F, ...
    penal, EminRatio);
caseNames = {'Baseline full-solid'; 'Optimized density'; 'Final clean redesign'};
shortNames = {'Baseline'; 'TO density'; 'Final clean'};
baselineCompliance = double(baselineData.compliance);
topologyCompliance = double(topOptData.compliance);
redesignCompliance = double(resultRedesign.compliance);
baselineMaxDisp = double(baselineData.maxDispMaterial);
topologyMaxDisp = double(topOptData.maxDispMaterial);
redesignMaxDisp = double(resultRedesign.maxDispMaterial);
baselineMaxVM = max(double(baselineData.vonMises(baselineData.materialElementIDs)));
topologyMaxVM = double(topOptData.maxVonMises);
redesignMaxVM = double(resultRedesign.maxVonMises);
caseVolume = [baselineVolume; topologyVolume; redesignedVolume];
caseCompliance = [baselineCompliance; topologyCompliance; redesignCompliance];
caseMaxDisp = [baselineMaxDisp; topologyMaxDisp; redesignMaxDisp];
caseMaxVM = [baselineMaxVM; topologyMaxVM; redesignMaxVM];
volumeRatioToBaseline = caseVolume / baselineVolume;
materialReductionPercent = 100 * (1 - volumeRatioToBaseline);
complianceRatioToBaseline = caseCompliance / baselineCompliance;
displacementRatioToBaseline = caseMaxDisp / baselineMaxDisp;
stiffnessRetentionPercent = 100 * baselineCompliance ./ caseCompliance;
comparisonTable = table(caseNames(:), caseVolume, volumeRatioToBaseline, materialReductionPercent, ...
    caseCompliance, complianceRatioToBaseline, stiffnessRetentionPercent, ...
    caseMaxDisp, displacementRatioToBaseline, caseMaxVM, ...
    'VariableNames', {'Case','MaterialVolume_mm3','VolumeRatioToBaseline','MaterialReduction_percent', ...
    'Compliance_Nmm','ComplianceRatioToBaseline','StiffnessRetention_percent', ...
    'MaxDisplacement_mm','DisplacementRatioToBaseline','MaxVonMises_MPa'});
redesignSummary = table;
redesignSummary.Item = {
    'Material case';
    'Young modulus E (MPa)';
    'Poisson ratio nu';
    'Thickness t (mm)';
    'Redesigned solid elements';
    'Redesigned solid nodes';
    'Fixed DOFs';
    'Free DOFs';
    'EminRatio';
    'penal';
    'Compliance C = F^T U (N mm)';
    'Total strain energy = 0.5 C (N mm)';
    'Max material displacement (mm)';
    'Node at max material displacement';
    'Max global displacement (mm)';
    'Max von Mises at element centre (MPa)';
    'Redesigned volume (mm^3)';
    'Volume ratio to baseline';
    'Material reduction vs baseline (%)';
    'Compliance ratio to baseline';
    'Displacement ratio to baseline';
    'Free DOF relative residual';
    'Total applied Fx (N)';
    'Total applied Fy (N)';
    'Total reaction Fx at fixed DOFs (N)';
    'Total reaction Fy at fixed DOFs (N)'};
redesignSummary.Value = {
    material.name;
    material.E_MPa;
    material.nu;
    material.thickness_mm;
    numel(solidElementIDs);
    numel(solidNodeIDs);
    numel(fixedDOFs);
    numel(freeDOFs);
    EminRatio;
    penal;
    resultRedesign.compliance;
    resultRedesign.strainEnergyTotal;
    resultRedesign.maxDispMaterial;
    resultRedesign.maxDispNode;
    resultRedesign.maxDispGlobal;
    resultRedesign.maxVonMises;
    redesignedVolume;
    redesignedVolume / baselineVolume;
    100 * (1 - redesignedVolume / baselineVolume);
    redesignCompliance / baselineCompliance;
    redesignMaxDisp / baselineMaxDisp;
    resultRedesign.freeDOFRelativeResidual;
    full(sum(F(1:2:end)));
    full(sum(F(2:2:end)));
    resultRedesign.reactionFx;
    resultRedesign.reactionFy};
writetable(comparisonTable, fullfile(outDir, 'validation_case_comparison.csv'));
writetable(redesignSummary, fullfile(outDir, 'validation_redesign_summary.csv'));
loadTable = struct2table(loadSpecsRedesign);
supportTable = struct2table(supportSpecsRedesign);
writetable(loadTable, fullfile(outDir, 'validation_loads_redesign.csv'));
writetable(supportTable, fullfile(outDir, 'validation_supports_redesign.csv'));
elementID = (1:nElements)';
elementIsRedesignSolid = xPhysRedesign >= solidThreshold;
elementStrainEnergy_Nmm = resultRedesign.elementStrainEnergy(:);
vonMises_MPa = resultRedesign.vonMises(:);
elementResults = table(elementID, elementIsRedesignSolid(:), xPhysRedesign(:), ...
    elementStrainEnergy_Nmm, vonMises_MPa);
writetable(elementResults, fullfile(outDir, 'validation_redesign_element_results.csv'));
fig = figure('Color','w', 'Position',[50 50 1600 900]);
imagesc(xCenters, yCenters, double(redesignMaskGrid)');
set(gca, 'YDir','normal'); axis equal tight; hold on; box on;
colormap(gca, gray(2));
cb = colorbar; ylabel(cb, '0 = void, 1 = solid');
plotBattery(batteryPoly, 'Battery cavity');
plotFrameMembers(P);
plotSupportsAndLoads(nodeCoords, supportSpecsRedesign, loadSpecsRedesign);
title('Final clean redesign validation setup with re-applied loads/supports');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig, fullfile(outDir, '27_validation_redesigned_setup.png'), exportResolution);
Ux = resultRedesign.U(1:2:end);
Uy = resultRedesign.U(2:2:end);
if resultRedesign.maxDispMaterial > eps && deformedShapeAutoScale
    modelSize = max([max(xEdges)-min(xEdges), max(yEdges)-min(yEdges)]);
    defScale = 0.10 * modelSize / resultRedesign.maxDispMaterial;
else
    defScale = manualDeformationScale;
end
deformedCoords = nodeCoords + defScale * [Ux, Uy];
fig = figure('Color','w', 'Position',[50 50 1600 900]);
hold on; axis equal; box on;
patch('Faces', elemNodes(solidElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceColor', [0.90 0.90 0.90], 'EdgeColor', [0.75 0.75 0.75], 'LineWidth', 0.2);
patch('Faces', elemNodes(solidElementIDs,:), 'Vertices', deformedCoords, ...
      'FaceVertexCData', resultRedesign.dispMag, 'FaceColor', 'interp', 'EdgeColor', edgeColor);
colorbar;
title(sprintf('Redesigned geometry deformed shape, scale = %.1f x', defScale));
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig, fullfile(outDir, '28_redesigned_deformed_shape.png'), exportResolution);
fig = figure('Color','w', 'Position',[50 50 1600 900]);
patch('Faces', elemNodes(solidElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', resultRedesign.dispMag, 'FaceColor', 'interp', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotSupportsAndLoads(nodeCoords, supportSpecsRedesign, loadSpecsRedesign);
title('Final clean redesign displacement magnitude');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig, fullfile(outDir, '29_redesigned_displacement_magnitude.png'), exportResolution);
fig = figure('Color','w', 'Position',[50 50 1600 900]);
patch('Faces', elemNodes(solidElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', resultRedesign.elementStrainEnergy(solidElementIDs), ...
      'FaceColor', 'flat', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotFrameMembers(P);
title('Final clean redesign element strain energy');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig, fullfile(outDir, '30_redesigned_element_strain_energy.png'), exportResolution);
fig = figure('Color','w', 'Position',[50 50 1600 900]);
patch('Faces', elemNodes(solidElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', resultRedesign.vonMises(solidElementIDs), ...
      'FaceColor', 'flat', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotFrameMembers(P);
title('Final clean redesign von Mises stress at Q4 element centre');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(fig, fullfile(outDir, '31_redesigned_von_mises_center.png'), exportResolution);
fig = figure('Color','w', 'Position',[50 50 1500 850]);
Y = [complianceRatioToBaseline, displacementRatioToBaseline, volumeRatioToBaseline];
bar(Y);
grid on; box on;
set(gca, 'XTickLabel', shortNames);
ylabel('Ratio to baseline full-solid case');
title('Final clean redesign validation: performance and material ratios');
legend({'Compliance ratio', 'Max displacement ratio', 'Material volume ratio'}, 'Location','northwest');
exportFigure(fig, fullfile(outDir, '32_validation_normalized_comparison.png'), exportResolution);
fig = figure('Color','w', 'Position',[50 50 1300 850]);
plot(volumeRatioToBaseline, complianceRatioToBaseline, 'ko-', 'LineWidth', 2, 'MarkerFaceColor','w', 'MarkerSize',8);
grid on; box on;
xlabel('Material volume ratio to baseline');
ylabel('Compliance ratio to baseline');
title('Compliance-volume trade-off');
for i = 1:numel(shortNames)
    text(volumeRatioToBaseline(i)+0.015, complianceRatioToBaseline(i), shortNames{i}, ...
        'FontWeight','bold', 'BackgroundColor','w');
end
exportFigure(fig, fullfile(outDir, '33_validation_compliance_volume_tradeoff.png'), exportResolution);
baselineMask = maskCode ~= 0;
topologyGrid = orientGrid(topOptData.xPhysGrid, nelx, nely, 'topOptData.xPhysGrid');
fig = figure('Color','w', 'Position',[50 50 1800 700]);
tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');
nexttile;
imagesc(xCenters, yCenters, double(baselineMask)'); set(gca,'YDir','normal'); axis equal tight; title('Baseline full-solid domain'); xlabel('x (mm)'); ylabel('y (mm)'); colormap(gca, gray(2));
nexttile;
imagesc(xCenters, yCenters, topologyGrid'); set(gca,'YDir','normal'); axis equal tight; title('Step 6 optimized density'); xlabel('x (mm)'); ylabel('y (mm)'); colormap(gca, flipud(gray(256))); colorbar;
nexttile;
imagesc(xCenters, yCenters, double(redesignMaskGrid)'); set(gca,'YDir','normal'); axis equal tight; title('Final clean redesign binary geometry'); xlabel('x (mm)'); ylabel('y (mm)'); colormap(gca, gray(2));
exportFigure(fig, fullfile(outDir, '34_validation_case_geometry_comparison.png'), exportResolution);
initialSkeletonMask = buildInitialSkeletonMask(xCenters, yCenters, P, batteryPoly, geometryData);
fig = figure('Color','w', 'Position',[50 50 1800 760]);
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact');
nexttile;
imagesc(xCenters, yCenters, double(initialSkeletonMask)');
set(gca,'YDir','normal'); axis equal tight; box on; colormap(gca, gray(2));
title('Initial reconstructed 2D CAD skeleton'); xlabel('x (mm)'); ylabel('y (mm)');
hold on; plotBattery(batteryPoly, 'Battery cavity'); plotFrameMembers(P);
nexttile;
imagesc(xCenters, yCenters, double(redesignMaskGrid)');
set(gca,'YDir','normal'); axis equal tight; box on; colormap(gca, gray(2));
title('Final clean 2D CAD redesign'); xlabel('x (mm)'); ylabel('y (mm)');
hold on; plotBattery(batteryPoly, 'Battery cavity'); plotFrameMembers(P);
sgtitle('Initial 2D CAD skeleton versus final clean redesign');
exportFigure(fig, fullfile(outDir, '35_final_clean_design_vs_initial_2d_cad.png'), exportResolution);

validationData = struct();
validationData.description = 'FEM validation of the final clean redesigned e-bike frame geometry';
validationData.settings = struct('solidThreshold', solidThreshold, ...
    'EminRatio', EminRatio, 'penal', penal, ...
    'reprojectLoadsAndSupportsToRedesign', reprojectLoadsAndSupportsToRedesign);
validationData.redesignedResult = resultRedesign;
validationData.comparisonTable = comparisonTable;
validationData.redesignSummary = redesignSummary;
validationData.loadSpecsRedesign = loadSpecsRedesign;
validationData.supportSpecsRedesign = supportSpecsRedesign;
validationData.solidElementIDs = solidElementIDs;
validationData.solidNodeIDs = solidNodeIDs;
validationData.xPhysRedesign = xPhysRedesign;
validationData.F = F;
validationData.fixedDOFs = fixedDOFs;
validationData.freeDOFs = freeDOFs;
save(fullfile(outDir, 'ebike_validation_step8.mat'), ...
    'validationData', 'redesignData', 'feMesh', 'geometryData', 'baselineData', 'topOptData', '-v7.3');
fprintf('\nFINAL CLEAN DESIGN - FEM validation completed successfully.\n');
fprintf('Output folder: %s\n\n', fullfile(pwd, outDir));
fprintf('Redesigned validation setup:\n');
fprintf('  Material case                     = %s\n', material.name);
fprintf('  E                                = %.3f MPa\n', material.E_MPa);
fprintf('  nu                               = %.3f\n', material.nu);
fprintf('  thickness                         = %.3f mm\n', material.thickness_mm);
fprintf('  Redesigned solid elements          = %d\n', numel(solidElementIDs));
fprintf('  Redesigned solid nodes             = %d\n', numel(solidNodeIDs));
fprintf('  Fixed DOFs                         = %d\n', numel(fixedDOFs));
fprintf('  Free DOFs                          = %d\n', numel(freeDOFs));
fprintf('  EminRatio                          = %.3e\n', EminRatio);
fprintf('  penal                              = %.3f\n\n', penal);
fprintf('Loads used in redesigned validation:\n');
for i = 1:numel(loadSpecsRedesign)
    p = nodeCoords(loadSpecsRedesign(i).nodeID, :);
    fprintf('  %-30s node %d at [%.2f, %.2f], Fx=%8.2f N, Fy=%8.2f N, snap distance=%.2f mm\n', ...
        loadSpecsRedesign(i).name, loadSpecsRedesign(i).nodeID, p(1), p(2), ...
        loadSpecsRedesign(i).Fx, loadSpecsRedesign(i).Fy, loadSpecsRedesign(i).snapDistance_mm);
end
fprintf('  Total applied Fx = %.3f N\n', full(sum(F(1:2:end))));
fprintf('  Total applied Fy = %.3f N\n\n', full(sum(F(2:2:end))));
fprintf('Supports used in redesigned validation:\n');
for i = 1:numel(supportSpecsRedesign)
    p = nodeCoords(supportSpecsRedesign(i).nodeID, :);
    fprintf('  %-24s node %d at [%.2f, %.2f], fixUx=%d, fixUy=%d, snap distance=%.2f mm\n', ...
        supportSpecsRedesign(i).name, supportSpecsRedesign(i).nodeID, p(1), p(2), ...
        supportSpecsRedesign(i).fixUx, supportSpecsRedesign(i).fixUy, supportSpecsRedesign(i).snapDistance_mm);
end
fprintf('  fixedDOFs = '); fprintf('%d ', fixedDOFs); fprintf('\n\n');
fprintf('Redesigned FEM results:\n');
fprintf('  Compliance C = F''*U              = %.6e N mm\n', resultRedesign.compliance);
fprintf('  Total strain energy = 0.5*C      = %.6e N mm\n', resultRedesign.strainEnergyTotal);
fprintf('  Max material displacement         = %.6e mm at node %d\n', resultRedesign.maxDispMaterial, resultRedesign.maxDispNode);
fprintf('  Max global displacement           = %.6e mm\n', resultRedesign.maxDispGlobal);
fprintf('  Max von Mises at element centre   = %.6e MPa\n', resultRedesign.maxVonMises);
fprintf('  Free DOF relative residual        = %.3e\n', resultRedesign.freeDOFRelativeResidual);
fprintf('  Reaction Fx at fixed DOFs         = %.6e N\n', resultRedesign.reactionFx);
fprintf('  Reaction Fy at fixed DOFs         = %.6e N\n\n', resultRedesign.reactionFy);
fprintf('Comparison to baseline full-solid domain:\n');
fprintf('  Redesign volume ratio             = %.4f\n', redesignedVolume / baselineVolume);
fprintf('  Material reduction                = %.2f %%\n', 100 * (1 - redesignedVolume / baselineVolume));
fprintf('  Compliance ratio                  = %.4f\n', redesignCompliance / baselineCompliance);
fprintf('  Max displacement ratio            = %.4f\n', redesignMaxDisp / baselineMaxDisp);
fprintf('  Stiffness retention estimate       = %.2f %%\n\n', 100 * baselineCompliance / redesignCompliance);
fprintf('Case comparison:\n');
for i = 1:numel(caseNames)
    fprintf('  %-28s | Vol ratio %.4f | C ratio %.4f | Disp ratio %.4f | Max VM %.2f MPa\n', ...
        caseNames{i}, volumeRatioToBaseline(i), complianceRatioToBaseline(i), ...
        displacementRatioToBaseline(i), caseMaxVM(i));
end
fprintf('\nFiles exported:\n');
fprintf('  27_validation_redesigned_setup.png\n');
fprintf('  28_redesigned_deformed_shape.png\n');
fprintf('  29_redesigned_displacement_magnitude.png\n');
fprintf('  30_redesigned_element_strain_energy.png\n');
fprintf('  31_redesigned_von_mises_center.png\n');
fprintf('  32_validation_normalized_comparison.png\n');
fprintf('  33_validation_compliance_volume_tradeoff.png\n');
fprintf('  34_validation_case_geometry_comparison.png\n');
fprintf('  35_final_clean_design_vs_initial_2d_cad.png\n');
fprintf('  validation_case_comparison.csv\n');
fprintf('  validation_redesign_summary.csv\n');
fprintf('  validation_loads_redesign.csv\n');
fprintf('  validation_supports_redesign.csv\n');
fprintf('  validation_redesign_element_results.csv\n');
fprintf('  ebike_validation_step8.mat\n\n');
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
function node = findClosestNode(nodeCoords, targetPoint, candidateNodeIDs)
    if nargin < 3 || isempty(candidateNodeIDs)
        candidateNodeIDs = (1:size(nodeCoords,1))';
    end
    candidateNodeIDs = candidateNodeIDs(:);
    candidateCoords = nodeCoords(candidateNodeIDs, :);
    d2 = (candidateCoords(:,1) - targetPoint(1)).^2 + ...
         (candidateCoords(:,2) - targetPoint(2)).^2;
    [~, idx] = min(d2);
    node = candidateNodeIDs(idx);
end
function result = solveStaticQ4Case(xPhys, solidElementIDs, solidNodeIDs, ...
        KE0, D, edofMat, elemNodes, nodeCoords, nDOF, fixedDOFs, freeDOFs, F, ...
        penal, EminRatio)
    nElements = size(edofMat, 1);
    xPhys = xPhys(:);
    stiffnessScale = EminRatio + (xPhys.^penal) * (1 - EminRatio);
    iK = reshape(kron(edofMat, ones(8,1))', 64*nElements, 1);
    jK = reshape(kron(edofMat, ones(1,8))', 64*nElements, 1);
    sK = reshape(KE0(:) * stiffnessScale(:)', 64*nElements, 1);
    K = sparse(iK, jK, sK, nDOF, nDOF);
    K = (K + K') / 2;
    U = zeros(nDOF, 1);
    U(freeDOFs) = K(freeDOFs, freeDOFs) \ F(freeDOFs);
    if any(~isfinite(U))
        error('FEM solution contains NaN or Inf. Check supports, loads, or redesigned connectivity.');
    end
    R = K * U - F;
    Ux = U(1:2:end);
    Uy = U(2:2:end);
    dispMag = sqrt(Ux.^2 + Uy.^2);
    [maxDispMaterial, idxLocal] = max(dispMag(solidNodeIDs));
    maxDispNode = solidNodeIDs(idxLocal);
    maxDispGlobal = max(dispMag);
    compliance = full(F' * U);
    strainEnergyTotal = 0.5 * compliance;
    Ue = U(edofMat);
    elementComplianceFullMaterial = sum((Ue * KE0) .* Ue, 2);
    elementComplianceActual = stiffnessScale .* elementComplianceFullMaterial;
    elementStrainEnergy = 0.5 * elementComplianceActual;
    n1 = elemNodes(1,1);
    n2 = elemNodes(1,2);
    n3 = elemNodes(1,3);
    n4 = elemNodes(1,4);
    coords4 = [nodeCoords(n1,:); nodeCoords(n2,:); nodeCoords(n3,:); nodeCoords(n4,:)];
    [~, Bcenter] = q4CenterBMatrix(coords4);
    stressCenter = zeros(nElements, 3);
    vonMises = zeros(nElements, 1);
    for e = solidElementIDs(:)'
        ue = Ue(e, :)';
        sigma = D * Bcenter * ue;
        stressCenter(e, :) = sigma';
        sx = sigma(1); sy = sigma(2); txy = sigma(3);
        vonMises(e) = sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2);
    end
    maxVonMises = max(vonMises(solidElementIDs));
    freeResidual = K(freeDOFs, freeDOFs) * U(freeDOFs) - F(freeDOFs);
    freeDOFRelativeResidual = norm(freeResidual) / max(1, norm(F(freeDOFs)));
    fixedOdd = fixedDOFs(mod(fixedDOFs,2) == 1);
    fixedEven = fixedDOFs(mod(fixedDOFs,2) == 0);
    reactionFx = full(sum(R(fixedOdd)));
    reactionFy = full(sum(R(fixedEven)));
    result = struct();
    result.U = U;
    result.R = R;
    result.K = K;
    result.stiffnessScale = stiffnessScale;
    result.dispMag = dispMag;
    result.compliance = compliance;
    result.strainEnergyTotal = strainEnergyTotal;
    result.maxDispMaterial = maxDispMaterial;
    result.maxDispNode = maxDispNode;
    result.maxDispGlobal = maxDispGlobal;
    result.elementComplianceFullMaterial = elementComplianceFullMaterial;
    result.elementComplianceActual = elementComplianceActual;
    result.elementStrainEnergy = elementStrainEnergy;
    result.stressCenter = stressCenter;
    result.vonMises = vonMises;
    result.maxVonMises = maxVonMises;
    result.freeDOFRelativeResidual = freeDOFRelativeResidual;
    result.reactionFx = reactionFx;
    result.reactionFy = reactionFy;
end
function [detJ, B] = q4CenterBMatrix(coords4)
    [detJ, B] = q4BMatrix(coords4, 0, 0);
end
function [detJ, B] = q4BMatrix(coords4, xi, eta)
    dN_dxi = 0.25 * [-(1-eta),  (1-eta),  (1+eta), -(1+eta)];
    dN_deta= 0.25 * [-(1-xi),  -(1+xi),   (1+xi),   (1-xi)];
    dN_nat = [dN_dxi; dN_deta];
    J = dN_nat * coords4;
    detJ = det(J);
    dN_global = J \ dN_nat;
    dNdx = dN_global(1,:);
    dNdy = dN_global(2,:);
    B = zeros(3,8);
    for a = 1:4
        col = 2*a - 1;
        B(1, col)     = dNdx(a);
        B(2, col + 1) = dNdy(a);
        B(3, col)     = dNdy(a);
        B(3, col + 1) = dNdx(a);
    end
end
function M = buildInitialSkeletonMask(xCenters, yCenters, P, batteryPoly, geometryData)
    [X, Y] = ndgrid(xCenters(:), yCenters(:));
    halfWidth = 17;
    if isfield(geometryData, 'referenceSkeleton') && isfield(geometryData.referenceSkeleton, 'tubeHalfWidth_mm')
        halfWidth = geometryData.referenceSkeleton.tubeHalfWidth_mm;
    end
    members = [rowPoint(P.RA); rowPoint(P.BB); rowPoint(P.RA); rowPoint(P.ST); ...
               rowPoint(P.BB); rowPoint(P.ST); rowPoint(P.ST); rowPoint(P.HT_top); ...
               rowPoint(P.BB); rowPoint(P.HT_bot); rowPoint(P.HT_bot); rowPoint(P.HT_top)];
    M = false(size(X));
    for ii = 1:2:size(members,1)
        M = M | segmentMaskLocal(X, Y, members(ii,:), members(ii+1,:), halfWidth);
    end
    if isfield(geometryData, 'anchor')
        a = geometryData.anchor;
        if isfield(a,'RA_radius'), M = M | pointMaskLocal(X,Y,rowPoint(P.RA),a.RA_radius); end
        if isfield(a,'BB_radius'), M = M | pointMaskLocal(X,Y,rowPoint(P.BB),a.BB_radius); end
        if isfield(a,'ST_radius'), M = M | pointMaskLocal(X,Y,rowPoint(P.ST),a.ST_radius); end
        if isfield(a,'HT_radius')
            M = M | pointMaskLocal(X,Y,rowPoint(P.HT_top),a.HT_radius);
            M = M | pointMaskLocal(X,Y,rowPoint(P.HT_bot),a.HT_radius);
        end
    end
    M = M & ~inpolygon(X, Y, batteryPoly(:,1), batteryPoly(:,2));
end
function mask = pointMaskLocal(X, Y, p, r)
    mask = (X-p(1)).^2 + (Y-p(2)).^2 <= r^2;
end
function mask = segmentMaskLocal(X, Y, p1, p2, halfWidth)
    vx = p2(1)-p1(1); vy = p2(2)-p1(2);
    wx = X-p1(1); wy = Y-p1(2);
    t = max(0, min(1, (wx*vx + wy*vy) ./ (vx^2 + vy^2)));
    projX = p1(1) + t*vx; projY = p1(2) + t*vy;
    mask = (X-projX).^2 + (Y-projY).^2 <= halfWidth^2;
end
function plotBattery(batteryPoly, displayName)
    patch('XData', batteryPoly(:,1), 'YData', batteryPoly(:,2), ...
          'FaceColor', [1 0.2 0.2], 'FaceAlpha', 0.12, ...
          'EdgeColor', [1 0 0], 'LineWidth', 2.0, 'DisplayName', displayName);
end
function plotFrameMembers(P)
    pairs = {P.RA, P.BB; P.RA, P.ST; P.BB, P.ST; P.ST, P.HT_top; P.BB, P.HT_bot; P.HT_bot, P.HT_top};
    for k = 1:size(pairs,1)
        p1 = rowPoint(pairs{k,1});
        p2 = rowPoint(pairs{k,2});
        plot([p1(1), p2(1)], [p1(2), p2(2)], 'k-', 'LineWidth', 1.2);
    end
end
function plotSupportsAndLoads(nodeCoords, supportSpecs, loadSpecs)
    for i = 1:numel(supportSpecs)
        nid = supportSpecs(i).nodeID;
        p = nodeCoords(nid,:);
        plot(p(1), p(2), 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        text(p(1) + 8, p(2) - 12, supportSpecs(i).name, ...
            'FontSize', 8, 'Interpreter', 'none', 'BackgroundColor','w');
    end
    for i = 1:numel(loadSpecs)
        nid = loadSpecs(i).nodeID;
        p = nodeCoords(nid,:);
        plot(p(1), p(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        drawLoadArrow(p, [loadSpecs(i).Fx, loadSpecs(i).Fy]);
        text(p(1) + 8, p(2) + 10, sprintf('%s\nFx=%.0f N, Fy=%.0f N', ...
            loadSpecs(i).name, loadSpecs(i).Fx, loadSpecs(i).Fy), ...
            'FontSize', 8, 'Interpreter', 'none', 'BackgroundColor','w');
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
function exportFigure(figHandle, filename, resolution)
    try
        exportgraphics(figHandle, filename, 'Resolution', resolution);
    catch
        saveas(figHandle, filename);
    end
end
