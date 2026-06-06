clear;
clc;
close all;

meshFile = fullfile(pwd, 'ebike_fe_mesh_outputs', 'ebike_fe_mesh_step3.mat');
outDir = fullfile(pwd, 'ebike_baseline_fem_outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

materialCase = 'paper_woven_equivalent';

customE_MPa = 70000;
customNu = 0.33;
customThickness_mm = 3.0;
analysisType = 'plane_stress';

penal = 3.0;
EminRatio = 1e-6;
doCleanupTinyActiveIslands = true;
minIslandElements = 20;
loads.rider_Fx_N   = 0;
loads.rider_Fy_N   = -850;
loads.battery_Fx_N = 0;
loads.battery_Fy_N = -50;
loads.motor_Fx_N   = 300;
loads.motor_Fy_N   = -300;
includeFrontImpact = false;
loads.frontImpact_Fx_N = -800;
loads.frontImpact_Fy_N = 0;
useRearAxleFixed = true;
useHeadTubeBottomRoller = true;
deformedShapeAutoScale = true;
manualDeformationScale = 30;
showElementEdgesInPlots = false;

if ~exist(meshFile, 'file')
    error(['Cannot find Step 3 mesh file:\n  %s\n\n', ...
           'Run generate_fe_mesh.m first. Make sure the folder ', ...
           'ebike_fe_mesh_outputs is inside your current MATLAB folder.'], meshFile);
end

load(meshFile, 'feMesh', 'bcData', 'geometryData');
P = geometryData.points;
batteryPoly = geometryData.batteryPoly;
nodeCoords = feMesh.nodeCoords;
elemNodes  = feMesh.elemNodes;
edofMat    = feMesh.elemDOFs;
maskCode   = feMesh.maskCode;
maskVector = maskCode(:);
nelx = feMesh.nelx;
nely = feMesh.nely;
nElements = feMesh.nElements;
nNodes = feMesh.nNodes;
nDOF = feMesh.nDOF;
xEdges = feMesh.xEdges;
yEdges = feMesh.yEdges;
xCenters = feMesh.xCenters;
yCenters = feMesh.yCenters;
XC = feMesh.XC;
YC = feMesh.YC;
dx = xEdges(2) - xEdges(1);
dy = yEdges(2) - yEdges(1);
material = selectMaterial(materialCase, customE_MPa, customNu, customThickness_mm);
E0 = material.E_MPa;
nu = material.nu;
thickness = material.thickness_mm;
if nu <= -1 || nu >= 0.5
    error('Poisson ratio must be between -1 and 0.5 for this linear elastic model.');
end
if E0 <= 0 || thickness <= 0
    error('Young modulus and thickness must be positive.');
end
D = planeStressD(E0, nu);
maskCodeClean = maskCode;
cleanupReport = struct('enabled', doCleanupTinyActiveIslands, ...
                       'removedActiveElements', 0, ...
                       'numMaterialComponentsBefore', NaN, ...
                       'numMaterialComponentsAfter', NaN);
if doCleanupTinyActiveIslands
    materialMaskBefore = maskCodeClean ~= 0;
    [compBefore, compSizesBefore] = connectedComponents4(materialMaskBefore);
    cleanupReport.numMaterialComponentsBefore = numel(compSizesBefore);
    for cid = 1:numel(compSizesBefore)
        compMask = (compBefore == cid);
        hasPassiveSolid = any(maskCodeClean(compMask) == 2);
        isTiny = compSizesBefore(cid) < minIslandElements;
        if isTiny && ~hasPassiveSolid
            tinyActive = compMask & (maskCodeClean == 1);
            cleanupReport.removedActiveElements = cleanupReport.removedActiveElements + nnz(tinyActive);
            maskCodeClean(tinyActive) = 0;
        end
    end
    materialMaskAfter = maskCodeClean ~= 0;
    [~, compSizesAfter] = connectedComponents4(materialMaskAfter);
    cleanupReport.numMaterialComponentsAfter = numel(compSizesAfter);
end
maskVectorClean = maskCodeClean(:);
activeElementIDs = find(maskVectorClean == 1);
passiveSolidElementIDs = find(maskVectorClean == 2);
passiveVoidElementIDs = find(maskVectorClean == 0);
materialElementIDs = find(maskVectorClean ~= 0);
materialNodeIDs = unique(elemNodes(materialElementIDs, :));
xPhys = zeros(nElements, 1);
xPhys(activeElementIDs) = 1.0;
xPhys(passiveSolidElementIDs) = 1.0;
stiffnessScale = EminRatio + (xPhys.^penal) * (1 - EminRatio);
supportSpecs = struct([]);
idx = 0;
if useRearAxleFixed
    idx = idx + 1;
    supportSpecs(idx).name = 'RA fixed support';
    supportSpecs(idx).targetPoint = P.RA;
    supportSpecs(idx).nodeID = findClosestNode(nodeCoords, P.RA, materialNodeIDs);
    supportSpecs(idx).fixUx = true;
    supportSpecs(idx).fixUy = true;
end
if useHeadTubeBottomRoller
    idx = idx + 1;
    supportSpecs(idx).name = 'HT_bot roller support';
    supportSpecs(idx).targetPoint = P.HT_bot;
    supportSpecs(idx).nodeID = findClosestNode(nodeCoords, P.HT_bot, materialNodeIDs);
    supportSpecs(idx).fixUx = false;
    supportSpecs(idx).fixUy = true;
end
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
batteryCenter = mean(batteryPoly, 1);
loadSpecs = struct([]);
loadSpecs(1).name = 'rider load at ST';
loadSpecs(1).targetPoint = P.ST;
loadSpecs(1).nodeID = findClosestNode(nodeCoords, P.ST, materialNodeIDs);
loadSpecs(1).Fx = loads.rider_Fx_N;
loadSpecs(1).Fy = loads.rider_Fy_N;
loadSpecs(2).name = 'battery load near down tube';
loadSpecs(2).targetPoint = batteryCenter;
loadSpecs(2).nodeID = findClosestNode(nodeCoords, batteryCenter, materialNodeIDs);
loadSpecs(2).Fx = loads.battery_Fx_N;
loadSpecs(2).Fy = loads.battery_Fy_N;
loadSpecs(3).name = 'motor/pedalling load at BB';
loadSpecs(3).targetPoint = P.BB;
loadSpecs(3).nodeID = findClosestNode(nodeCoords, P.BB, materialNodeIDs);
loadSpecs(3).Fx = loads.motor_Fx_N;
loadSpecs(3).Fy = loads.motor_Fy_N;
if includeFrontImpact
    loadSpecs(4).name = 'optional front impact at HT_top';
    loadSpecs(4).targetPoint = P.HT_top;
    loadSpecs(4).nodeID = findClosestNode(nodeCoords, P.HT_top, materialNodeIDs);
    loadSpecs(4).Fx = loads.frontImpact_Fx_N;
    loadSpecs(4).Fy = loads.frontImpact_Fy_N;
end
F = sparse(nDOF, 1);
for i = 1:numel(loadSpecs)
    nid = loadSpecs(i).nodeID;
    F(2*nid - 1) = F(2*nid - 1) + loadSpecs(i).Fx;
    F(2*nid)     = F(2*nid)     + loadSpecs(i).Fy;
end
coordsRect = [0, 0;
              dx, 0;
              dx, dy;
              0, dy];
KE0 = q4ElementStiffness(coordsRect, D, thickness);
KEsymError = norm(KE0 - KE0', 'fro') / max(1, norm(KE0, 'fro'));
iK = reshape(kron(edofMat, ones(8,1))', 64*nElements, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nElements, 1);
sK = reshape(KE0(:) * stiffnessScale(:)', 64*nElements, 1);
K = sparse(iK, jK, sK, nDOF, nDOF);
K = (K + K') / 2;
Kff = K(freeDOFs, freeDOFs);
Ff = F(freeDOFs);
U = zeros(nDOF, 1);
U(freeDOFs) = Kff \ Ff;
if any(~isfinite(U))
    error('The FEM solution contains NaN or Inf. Check supports, loads, and EminRatio.');
end
R = K * U - F;
Ux = U(1:2:end);
Uy = U(2:2:end);
dispMag = sqrt(Ux.^2 + Uy.^2);
maxDispGlobal = max(dispMag);
[maxDispMaterial, idxMatMaxLocal] = max(dispMag(materialNodeIDs));
maxDispNode = materialNodeIDs(idxMatMaxLocal);
compliance = full(F' * U);
externalWork = compliance;
strainEnergyTotal = 0.5 * compliance;
Ue = U(edofMat);
elementComplianceFullMaterial = sum((Ue * KE0) .* Ue, 2);
elementComplianceActual = stiffnessScale .* elementComplianceFullMaterial;
elementStrainEnergy = 0.5 * elementComplianceActual;
[~, Bcenter] = q4CenterBMatrix(coordsRect);
stressCenter = zeros(nElements, 3);
vonMises = zeros(nElements, 1);
for e = materialElementIDs(:)'
    ue = Ue(e, :)';
    sigma = D * Bcenter * ue;
    stressCenter(e, :) = sigma';
    sx = sigma(1); sy = sigma(2); txy = sigma(3);
    vonMises(e) = sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2);
end
maxVonMises = max(vonMises(materialElementIDs));
reactionAtFixedDOFs = full(R(fixedDOFs));
elementArea = dx * dy;
elementVolume = elementArea * thickness;
materialAreaBaseline = numel(materialElementIDs) * elementArea;
materialVolumeBaseline = numel(materialElementIDs) * elementVolume;
designAreaActivePlusSolid = nnz(maskVectorClean ~= 0) * elementArea;
activeArea = numel(activeElementIDs) * elementArea;
passiveSolidArea = numel(passiveSolidElementIDs) * elementArea;
voidArea = numel(passiveVoidElementIDs) * elementArea;
edgeColor = 'none';
if showElementEdgesInPlots
    edgeColor = [0.70 0.70 0.70];
end
figA = figure('Color','w','Name','Step 4 loads and supports');
hold on; axis equal; box on;
imagesc(xCenters, yCenters, maskCodeClean);
set(gca, 'YDir', 'normal');
colormap(gca, [0.96 0.96 0.96; 0.75 0.85 1.00; 0.20 0.45 0.95]);
cb = colorbar;
cb.Ticks = [0, 1, 2];
cb.TickLabels = {'void/outside','active solid baseline','passive solid'};
title('Step 4: material, loads and supports for baseline FEM');
xlabel('x position (mm)'); ylabel('y position (mm)');
plotFrameMembers(P);
plot([batteryPoly(:,1); batteryPoly(1,1)], [batteryPoly(:,2); batteryPoly(1,2)], 'r-', 'LineWidth', 2);
plotSupportsAndLoads(nodeCoords, supportSpecs, loadSpecs);
exportFigure(figA, fullfile(outDir, '07_step4_material_loads_supports.png'));
if maxDispMaterial > eps && deformedShapeAutoScale
    modelSize = max([max(xEdges)-min(xEdges), max(yEdges)-min(yEdges)]);
    defScale = 0.10 * modelSize / maxDispMaterial;
else
    defScale = manualDeformationScale;
end
deformedCoords = nodeCoords + defScale * [Ux, Uy];
figB = figure('Color','w','Name','Baseline deformed shape');
hold on; axis equal; box on;
patch('Faces', elemNodes(materialElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceColor', [0.90 0.90 0.90], 'EdgeColor', [0.75 0.75 0.75], 'LineWidth', 0.2);
patch('Faces', elemNodes(materialElementIDs,:), 'Vertices', deformedCoords, ...
      'FaceVertexCData', dispMag, 'FaceColor', 'interp', 'EdgeColor', edgeColor);
colorbar;
title(sprintf('Baseline FEM deformed shape, scale = %.1f x', defScale));
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(figB, fullfile(outDir, '08_baseline_deformed_shape.png'));
figC = figure('Color','w','Name','Baseline displacement magnitude');
patch('Faces', elemNodes(materialElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', dispMag, 'FaceColor', 'interp', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotSupportsAndLoads(nodeCoords, supportSpecs, loadSpecs);
title('Baseline FEM displacement magnitude on material domain');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(figC, fullfile(outDir, '09_baseline_displacement_magnitude.png'));
figD = figure('Color','w','Name','Baseline element strain energy');
patch('Faces', elemNodes(materialElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', elementStrainEnergy(materialElementIDs), ...
      'FaceColor', 'flat', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotFrameMembers(P);
title('Baseline FEM element strain energy');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(figD, fullfile(outDir, '10_baseline_element_strain_energy.png'));
figE = figure('Color','w','Name','Baseline von Mises stress');
patch('Faces', elemNodes(materialElementIDs,:), 'Vertices', nodeCoords, ...
      'FaceVertexCData', vonMises(materialElementIDs), ...
      'FaceColor', 'flat', 'EdgeColor', edgeColor);
hold on; axis equal tight; box on; colorbar;
plotFrameMembers(P);
title('Baseline FEM von Mises stress at Q4 element centre');
xlabel('x position (mm)'); ylabel('y position (mm)');
exportFigure(figE, fullfile(outDir, '11_baseline_von_mises_center.png'));
summary = table;
summary.Item = {
    'Material case';
    'Young modulus E (MPa)';
    'Poisson ratio nu';
    'Thickness t (mm)';
    'Analysis type';
    'Mesh nelx';
    'Mesh nely';
    'Total elements';
    'Active design elements';
    'Passive solid elements';
    'Passive void/outside elements';
    'Removed tiny active islands';
    'Total DOFs';
    'Fixed DOFs';
    'Free DOFs';
    'Compliance C = F^T U (N mm)';
    'Total strain energy = 0.5 C (N mm)';
    'Max material displacement (mm)';
    'Node at max material displacement';
    'Max global displacement including void nodes (mm)';
    'Max von Mises at element centre (MPa)';
    'Baseline material area (mm^2)';
    'Baseline material volume (mm^3)';
    'Element stiffness symmetry error'};
summary.Value = {
    material.name;
    E0;
    nu;
    thickness;
    analysisType;
    nelx;
    nely;
    nElements;
    numel(activeElementIDs);
    numel(passiveSolidElementIDs);
    numel(passiveVoidElementIDs);
    cleanupReport.removedActiveElements;
    nDOF;
    numel(fixedDOFs);
    numel(freeDOFs);
    compliance;
    strainEnergyTotal;
    maxDispMaterial;
    maxDispNode;
    maxDispGlobal;
    maxVonMises;
    materialAreaBaseline;
    materialVolumeBaseline;
    KEsymError};
try
    writetable(summary, fullfile(outDir, 'baseline_fem_summary.csv'));
catch
    warning('Could not write baseline_fem_summary.csv using writetable.');
end
loadTable = struct2table(loadSpecs);
supportTable = struct2table(supportSpecs);
try
    writetable(loadTable, fullfile(outDir, 'step4_loads.csv'));
    writetable(supportTable, fullfile(outDir, 'step4_supports.csv'));
catch
    warning('Could not write load/support CSV tables.');
end
baselineData = struct();
baselineData.material = material;
baselineData.loads = loads;
baselineData.loadSpecs = loadSpecs;
baselineData.supportSpecs = supportSpecs;
baselineData.F = F;
baselineData.fixedDOFs = fixedDOFs;
baselineData.freeDOFs = freeDOFs;
baselineData.U = U;
baselineData.R = R;
baselineData.compliance = compliance;
baselineData.strainEnergyTotal = strainEnergyTotal;
baselineData.maxDispMaterial = maxDispMaterial;
baselineData.maxDispNode = maxDispNode;
baselineData.maxDispGlobal = maxDispGlobal;
baselineData.elementStrainEnergy = elementStrainEnergy;
baselineData.elementComplianceActual = elementComplianceActual;
baselineData.stressCenter = stressCenter;
baselineData.vonMises = vonMises;
baselineData.xPhys = xPhys;
baselineData.stiffnessScale = stiffnessScale;
baselineData.KE0 = KE0;
baselineData.D = D;
baselineData.penal = penal;
baselineData.EminRatio = EminRatio;
baselineData.maskCodeClean = maskCodeClean;
baselineData.maskVectorClean = maskVectorClean;
baselineData.activeElementIDs = activeElementIDs;
baselineData.passiveSolidElementIDs = passiveSolidElementIDs;
baselineData.passiveVoidElementIDs = passiveVoidElementIDs;
baselineData.materialElementIDs = materialElementIDs;
baselineData.materialNodeIDs = materialNodeIDs;
baselineData.cleanupReport = cleanupReport;
baselineData.elementArea = elementArea;
baselineData.elementVolume = elementVolume;
baselineData.materialAreaBaseline = materialAreaBaseline;
baselineData.materialVolumeBaseline = materialVolumeBaseline;
baselineData.reactionAtFixedDOFs = reactionAtFixedDOFs;
save(fullfile(outDir, 'ebike_baseline_fem_step4_5.mat'), ...
     'feMesh', 'geometryData', 'baselineData', 'summary', 'loadTable', 'supportTable');
fprintf('\nSTEP 4 + STEP 5 completed successfully.\n');
fprintf('Output folder: %s\n', outDir);
fprintf('\nMaterial model:\n');
fprintf('  Case       = %s\n', material.name);
fprintf('  E          = %.3f MPa\n', E0);
fprintf('  nu         = %.3f\n', nu);
fprintf('  thickness  = %.3f mm\n', thickness);
fprintf('  Note       = %s\n', material.note);
fprintf('\nMesh and mask after cleanup:\n');
fprintf('  nelx x nely                  = %d x %d\n', nelx, nely);
fprintf('  Total elements               = %d\n', nElements);
fprintf('  Active baseline solid elems  = %d\n', numel(activeElementIDs));
fprintf('  Passive solid anchor elems   = %d\n', numel(passiveSolidElementIDs));
fprintf('  Passive void/outside elems   = %d\n', numel(passiveVoidElementIDs));
fprintf('  Tiny active elements removed = %d\n', cleanupReport.removedActiveElements);
fprintf('\nSupports:\n');
for i = 1:numel(supportSpecs)
    p = nodeCoords(supportSpecs(i).nodeID, :);
    fprintf('  %-24s node %d at [%.2f, %.2f], fixUx=%d, fixUy=%d\n', ...
        supportSpecs(i).name, supportSpecs(i).nodeID, p(1), p(2), ...
        supportSpecs(i).fixUx, supportSpecs(i).fixUy);
end
fprintf('  fixedDOFs = '); fprintf('%d ', fixedDOFs); fprintf('\n');
fprintf('\nLoads:\n');
for i = 1:numel(loadSpecs)
    p = nodeCoords(loadSpecs(i).nodeID, :);
    fprintf('  %-30s node %d at [%.2f, %.2f], Fx=%.1f N, Fy=%.1f N\n', ...
        loadSpecs(i).name, loadSpecs(i).nodeID, p(1), p(2), ...
        loadSpecs(i).Fx, loadSpecs(i).Fy);
end
fprintf('  Total applied Fx = %.3f N\n', full(sum(F(1:2:end))));
fprintf('  Total applied Fy = %.3f N\n', full(sum(F(2:2:end))));
fprintf('\nBaseline FEM results:\n');
fprintf('  Compliance C = F''*U              = %.6e N mm\n', compliance);
fprintf('  Total strain energy = 0.5*C      = %.6e N mm\n', strainEnergyTotal);
fprintf('  Max material displacement         = %.6e mm at node %d\n', maxDispMaterial, maxDispNode);
fprintf('  Max global displacement           = %.6e mm\n', maxDispGlobal);
fprintf('  Max von Mises at element centre   = %.6e MPa\n', maxVonMises);
fprintf('  Element stiffness symmetry error  = %.3e\n', KEsymError);
fprintf('\nFiles exported:\n');
fprintf('  07_step4_material_loads_supports.png\n');
fprintf('  08_baseline_deformed_shape.png\n');
fprintf('  09_baseline_displacement_magnitude.png\n');
fprintf('  10_baseline_element_strain_energy.png\n');
fprintf('  11_baseline_von_mises_center.png\n');
fprintf('  baseline_fem_summary.csv\n');
fprintf('  step4_loads.csv\n');
fprintf('  step4_supports.csv\n');
fprintf('  ebike_baseline_fem_step4_5.mat\n\n');
function material = selectMaterial(materialCase, customE, customNu, customT)
    switch lower(materialCase)
        case 'paper_woven_equivalent'
            material.name = 'paper_woven_equivalent';
            material.E_MPa = 61340;
            material.nu = 0.30;
            material.thickness_mm = 1.5;
            material.note = ['Simplified isotropic equivalent using Yasin et al. woven prepreg ', ...
                             'longitudinal E1 and v12. Not a full laminate model.'];
        case 'paper_ud_equivalent'
            material.name = 'paper_ud_equivalent';
            material.E_MPa = 121000;
            material.nu = 0.27;
            material.thickness_mm = 1.5;
            material.note = ['Simplified isotropic equivalent using Yasin et al. UD prepreg ', ...
                             'longitudinal E1 and v12. Not a full laminate model.'];
        case 'aluminium6061'
            material.name = 'aluminium6061';
            material.E_MPa = 70000;
            material.nu = 0.33;
            material.thickness_mm = 3.0;
            material.note = 'Common isotropic aluminium approximation; easier, but not the Yasin paper material.';
        case 'custom'
            material.name = 'custom';
            material.E_MPa = customE;
            material.nu = customNu;
            material.thickness_mm = customT;
            material.note = 'User-defined isotropic material.';
        otherwise
            error('Unknown materialCase. Use paper_woven_equivalent, paper_ud_equivalent, aluminium6061, or custom.');
    end
end
function D = planeStressD(E, nu)
    D = E/(1 - nu^2) * [1,  nu, 0;
                        nu, 1,  0;
                        0,  0,  (1 - nu)/2];
end
function KE = q4ElementStiffness(coords4, D, thickness)
    gp = [-1/sqrt(3), 1/sqrt(3)];
    KE = zeros(8,8);
    for i = 1:2
        xi = gp(i);
        for j = 1:2
            eta = gp(j);
            [detJ, B] = q4BMatrix(coords4, xi, eta);
            if detJ <= 0
                error('Q4 element has non-positive Jacobian determinant. Check element node ordering.');
            end
            KE = KE + (B' * D * B) * detJ * thickness;
        end
    end
    KE = (KE + KE') / 2;
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
function node = findClosestNode(nodeCoords, targetPoint, candidateNodeIDs)
    if nargin < 3 || isempty(candidateNodeIDs)
        candidateNodeIDs = (1:size(nodeCoords,1))';
    end
    candidateCoords = nodeCoords(candidateNodeIDs, :);
    d2 = (candidateCoords(:,1) - targetPoint(1)).^2 + ...
         (candidateCoords(:,2) - targetPoint(2)).^2;
    [~, idx] = min(d2);
    node = candidateNodeIDs(idx);
end
function [labels, sizes] = connectedComponents4(binaryMask)
    [nr, nc] = size(binaryMask);
    labels = zeros(nr, nc);
    sizes = [];
    comp = 0;
    maxQueue = nr * nc;
    qR = zeros(maxQueue, 1);
    qC = zeros(maxQueue, 1);
    for r = 1:nr
        for c = 1:nc
            if binaryMask(r,c) && labels(r,c) == 0
                comp = comp + 1;
                head = 1;
                tail = 1;
                qR(tail) = r;
                qC(tail) = c;
                labels(r,c) = comp;
                count = 0;
                while head <= tail
                    rr = qR(head);
                    cc = qC(head);
                    head = head + 1;
                    count = count + 1;
                    neigh = [rr-1, cc;
                             rr+1, cc;
                             rr, cc-1;
                             rr, cc+1];
                    for k = 1:4
                        r2 = neigh(k,1);
                        c2 = neigh(k,2);
                        if r2 >= 1 && r2 <= nr && c2 >= 1 && c2 <= nc
                            if binaryMask(r2,c2) && labels(r2,c2) == 0
                                tail = tail + 1;
                                qR(tail) = r2;
                                qC(tail) = c2;
                                labels(r2,c2) = comp;
                            end
                        end
                    end
                end
                sizes(comp,1) = count;
            end
        end
    end
end
function plotFrameMembers(P)
    pairs = {P.RA, P.BB; P.RA, P.ST; P.BB, P.ST; P.ST, P.HT_top; P.BB, P.HT_bot; P.HT_bot, P.HT_top};
    for k = 1:size(pairs,1)
        p1 = pairs{k,1};
        p2 = pairs{k,2};
        plot([p1(1), p2(1)], [p1(2), p2(2)], 'k-', 'LineWidth', 1.5);
    end
end
function plotSupportsAndLoads(nodeCoords, supportSpecs, loadSpecs)
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
