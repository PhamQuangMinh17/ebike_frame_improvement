clear;
clc;
close all;

wheelbase       = 1086.1;
headTubeAngle   = 69;
seatTubeAngle   = 73;
seatTubeLength  = 500;
seatStayLength  = 446.6;
topTubeLength   = 607.6;
chainStayLength = 430;
headTubeLength  = 105;
plotTubeWidth = 18;
referenceTubeWidth_mm = 34;
referenceTubeHalfWidth_mm = referenceTubeWidth_mm/2;
nelx = 160;
nely = 90;
margin = 80;
anchor.RA_radius = 45;
anchor.BB_radius = 60;
anchor.ST_radius = 42;
anchor.HT_radius= 45;
anchor.HeadTube_halfWidth = 28;
battery.length_mm = 330;
battery.width_mm  = 65;
battery.centerFractionAlongDownTube = 0.55;
battery.offsetTowardInside_mm = 35;
loadInfo.riderLoad_N  = -800;
loadInfo.batteryLoad_N = -80;
loadInfo.motorLoad_N = -150;

outDir = fullfile(pwd, 'ebike_geometry_outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
RA = [0, 0];
FA = [wheelbase, 0];
seatTheta = deg2rad(seatTubeAngle);
headTheta = deg2rad(headTubeAngle);
objectiveAlpha = @(alpha) norm(chainStayLength*[cos(alpha), sin(alpha)] + ...
    seatTubeLength*[-cos(seatTheta), sin(seatTheta)] - RA) - seatStayLength;
try
    alpha = fzero(objectiveAlpha, [-pi/2, 0]);
catch
    alpha = deg2rad(-16.2);
end

BB = chainStayLength * [cos(alpha), sin(alpha)];
ST = BB + seatTubeLength * [-cos(seatTheta), sin(seatTheta)];
HT_top = ST + [topTubeLength, 0];
HT_bot = HT_top + headTubeLength * [cos(headTheta), -sin(headTheta)];
P.RA = RA;
P.BB = BB;
P.ST = ST;
P.HT_top = HT_top;
P.HT_bot = HT_bot;
P.FA = FA;

memberNames = {'chain stay','seat stay','seat tube','top tube','down tube','head tube'};
members = [
    RA;      BB;
    RA;      ST;
    BB;      ST;
    ST;      HT_top;
    BB;      HT_bot;
    HT_bot;  HT_top
];
rearTri = [RA; BB; ST];
mainQuad = [BB; ST; HT_top; HT_bot];
downTubeStart = BB;
downTubeEnd   = HT_bot;
batteryPoly = makeOrientedRectangle(downTubeStart, downTubeEnd, ...
    battery.centerFractionAlongDownTube, battery.length_mm, battery.width_mm, ...
    battery.offsetTowardInside_mm);
allPts = [RA; BB; ST; HT_top; HT_bot; FA];
xmin = min(allPts(:,1)) - margin;
xmax = max(allPts(:,1)) + margin;
ymin = min(allPts(:,2)) - margin;
ymax = max(allPts(:,2)) + margin;
xEdges = linspace(xmin, xmax, nelx+1);
yEdges = linspace(ymin, ymax, nely+1);
xCenters = 0.5*(xEdges(1:end-1) + xEdges(2:end));
yCenters = 0.5*(yEdges(1:end-1) + yEdges(2:end));
[XC, YC] = meshgrid(xCenters, yCenters);
insideRear = inpolygon(XC, YC, rearTri(:,1), rearTri(:,2));
insideMain = inpolygon(XC, YC, mainQuad(:,1), mainQuad(:,2));
activeMask = insideRear | insideMain;
insideBattery = inpolygon(XC, YC, batteryPoly(:,1), batteryPoly(:,2));
batteryVoidMask = activeMask & insideBattery;

passiveSolidMask = false(size(activeMask));
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, RA, anchor.RA_radius);
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, BB, anchor.BB_radius);
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, ST, anchor.ST_radius);
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, HT_top, anchor.HT_radius);
passiveSolidMask = passiveSolidMask | pointRadiusMask(XC, YC, HT_bot, anchor.HT_radius);
passiveSolidMask = passiveSolidMask | segmentDistanceMask(XC, YC, HT_bot, HT_top, anchor.HeadTube_halfWidth);
passiveSolidMask = activeMask & passiveSolidMask & ~batteryVoidMask;
outsideDesignMask = ~activeMask;
passiveVoidMask = outsideDesignMask | batteryVoidMask;

maskCode = zeros(size(activeMask));
maskCode(activeMask) = 1;
maskCode(passiveVoidMask) = 0;
maskCode(passiveSolidMask) = 2;
referenceSkeletonMask = false(size(activeMask));
for i = 1:2:size(members,1)
    referenceSkeletonMask = referenceSkeletonMask | segmentDistanceMask(XC, YC, members(i,:), members(i+1,:), referenceTubeHalfWidth_mm);
end
referenceSkeletonMask = (referenceSkeletonMask | passiveSolidMask) & activeMask & ~batteryVoidMask;
fig1 = figure('Color','w','Name','Baseline e-bike frame skeleton');
hold on; axis equal; grid on; box on;
title('Baseline e-bike frame skeleton from Yasin et al. (2023) dimensions');
xlabel('x position (mm)'); ylabel('y position (mm)');
for i = 1:2:size(members,1)
    p1 = members(i,:); p2 = members(i+1,:);
    plot([p1(1),p2(1)], [p1(2),p2(2)], 'k-', 'LineWidth', 3);
end
plotPointWithLabel(RA, 'RA');
plotPointWithLabel(BB, 'BB / motor');
plotPointWithLabel(ST, 'ST');
plotPointWithLabel(HT_top, 'HT top');
plotPointWithLabel(HT_bot, 'HT bot');
plotPointWithLabel(FA, 'FA ref');
text(xmin+20, ymax-25, sprintf(['Dimension checks (mm):\n', ...
    'RA-BB chain stay = %.1f\n', ...
    'RA-ST seat stay = %.1f\n', ...
    'BB-ST seat tube = %.1f\n', ...
    'ST-HT top tube = %.1f\n', ...
    'HT tube = %.1f\n', ...
    'Wheelbase = %.1f'], ...
    norm(BB-RA), norm(ST-RA), norm(ST-BB), norm(HT_top-ST), norm(HT_top-HT_bot), wheelbase), ...
    'VerticalAlignment','top', 'FontName','Consolas');
xlim([xmin xmax]); ylim([ymin ymax]);
exportFigure(fig1, fullfile(outDir, '01_keypoints_baseline_skeleton.png'));
figRef = figure('Color','w','Name','Reference 2D skeleton CAD');
imagesc(xCenters, yCenters, double(referenceSkeletonMask));
set(gca,'YDir','normal'); axis equal tight; grid on; box on;
colormap(gca, gray(2));
title(sprintf('Reference 2D frame skeleton CAD, equivalent member width = %.0f mm', referenceTubeWidth_mm));
xlabel('x position (mm)'); ylabel('y position (mm)');
hold on;
for i = 1:2:size(members,1)
    p1 = members(i,:); p2 = members(i+1,:);
    plot([p1(1),p2(1)], [p1(2),p2(2)], 'k-', 'LineWidth', 1.2);
end
patch(batteryPoly(:,1), batteryPoly(:,2), [1 1 1], 'FaceAlpha', 0.20, 'EdgeColor', 'r', 'LineWidth', 2);
text(mean(batteryPoly(:,1)), mean(batteryPoly(:,2)), {'battery','cavity'}, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','none', 'Color','r');
plotPointWithLabel(RA, 'RA'); plotPointWithLabel(BB, 'BB'); plotPointWithLabel(ST, 'ST');
plotPointWithLabel(HT_top, 'HT top'); plotPointWithLabel(HT_bot, 'HT bot');
exportFigure(figRef, fullfile(outDir, '01b_reference_skeleton_2d_cad.png'));
fig2 = figure('Color','w','Name','Active design region and passive regions');
hold on; axis equal; grid on; box on;
title('MATLAB-generated CAD / geometric model for topology optimization');
xlabel('x position (mm)'); ylabel('y position (mm)');
patch(rearTri(:,1), rearTri(:,2), [0.85 0.85 0.85], 'FaceAlpha', 0.45, 'EdgeColor', 'k');
patch(mainQuad(:,1), mainQuad(:,2), [0.85 0.85 0.85], 'FaceAlpha', 0.45, 'EdgeColor', 'k');
patch(batteryPoly(:,1), batteryPoly(:,2), [1 1 1], 'FaceAlpha', 0.95, 'EdgeColor', 'r', 'LineWidth', 2);
text(mean(batteryPoly(:,1)), mean(batteryPoly(:,2)), {'battery cavity','passive void'}, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','none');
drawCircle(RA, anchor.RA_radius, 'RA anchor');
drawCircle(BB, anchor.BB_radius, 'BB / motor anchor');
drawCircle(ST, anchor.ST_radius, 'ST anchor');
drawCircle(HT_top, anchor.HT_radius, 'HT top anchor');
drawCircle(HT_bot, anchor.HT_radius, 'HT bot anchor');
for i = 1:2:size(members,1)
    p1 = members(i,:); p2 = members(i+1,:);
    plot([p1(1),p2(1)], [p1(2),p2(2)], 'k-', 'LineWidth', 2);
end
legend({'Rear active region','Main active region','Battery passive void'}, 'Location','southoutside');
xlim([xmin xmax]); ylim([ymin ymax]);
exportFigure(fig2, fullfile(outDir, '02_design_domain_passive_regions.png'));
fig3 = figure('Color','w','Name','Element mask for topology optimization');
imagesc(xCenters, yCenters, maskCode);
set(gca,'YDir','normal'); axis equal tight; grid on; box on;
title('Element mask: 0 = passive void, 1 = active design, 2 = passive solid');
xlabel('x position (mm)'); ylabel('y position (mm)');
colorbar;
hold on;
plot([RA(1),BB(1),ST(1),HT_top(1),HT_bot(1),BB(1)], ...
     [RA(2),BB(2),ST(2),HT_top(2),HT_bot(2),BB(2)], 'k-', 'LineWidth', 1.5);
plot(batteryPoly(:,1), batteryPoly(:,2), 'r-', 'LineWidth', 2);
exportFigure(fig3, fullfile(outDir, '03_grid_masks_for_topology_optimization.png'));
geometryData = struct();
geometryData.reference = 'Yasin et al. (2023) e-bike frame dimensions; reconstructed 2D MATLAB model';
geometryData.units = 'mm, N, MPa';
geometryData.dimensions = struct('wheelbase',wheelbase, 'headTubeAngle',headTubeAngle, ...
    'seatTubeAngle',seatTubeAngle, 'seatTubeLength',seatTubeLength, ...
    'seatStayLength',seatStayLength, 'topTubeLength',topTubeLength, ...
    'chainStayLength',chainStayLength, 'headTubeLength',headTubeLength);
geometryData.points = P;
geometryData.memberNames = memberNames;
geometryData.members = members;
geometryData.rearTri = rearTri;
geometryData.mainQuad = mainQuad;
geometryData.batteryPoly = batteryPoly;
geometryData.grid = struct('nelx',nelx,'nely',nely,'xEdges',xEdges,'yEdges',yEdges, ...
    'xCenters',xCenters,'yCenters',yCenters,'XC',XC,'YC',YC);
geometryData.masks = struct('activeMask',activeMask, 'passiveVoidMask',passiveVoidMask, ...
    'batteryVoidMask',batteryVoidMask, 'passiveSolidMask',passiveSolidMask, 'maskCode',maskCode);
geometryData.referenceSkeleton = struct('mask',referenceSkeletonMask, ...
    'tubeWidth_mm',referenceTubeWidth_mm, 'tubeHalfWidth_mm',referenceTubeHalfWidth_mm, ...
    'description','Reconstructed 2D member-based skeleton used to visualise the initial frame CAD before the filled topology design domain.');
geometryData.loadInfo = loadInfo;
geometryData.anchor = anchor;
geometryData.battery = battery;
save(fullfile(outDir, 'ebike_geometry_yasin2023.mat'), 'geometryData');
fprintf('\nMATLAB CAD / geometric model generated successfully.\n');
fprintf('Output folder: %s\n', outDir);
fprintf('Key reconstructed points [x, y] in mm:\n');
fprintf('  RA     = [%8.2f, %8.2f]\n', RA(1), RA(2));
fprintf('  BB     = [%8.2f, %8.2f]\n', BB(1), BB(2));
fprintf('  ST     = [%8.2f, %8.2f]\n', ST(1), ST(2));
fprintf('  HT_top = [%8.2f, %8.2f]\n', HT_top(1), HT_top(2));
fprintf('  HT_bot = [%8.2f, %8.2f]\n', HT_bot(1), HT_bot(2));
fprintf('  FA     = [%8.2f, %8.2f]\n', FA(1), FA(2));
fprintf('\nFiles exported:\n');
fprintf('  01_keypoints_baseline_skeleton.png\n');
fprintf('  01b_reference_skeleton_2d_cad.png\n');
fprintf('  02_design_domain_passive_regions.png\n');
fprintf('  03_grid_masks_for_topology_optimization.png\n');
fprintf('  ebike_geometry_yasin2023.mat\n\n');
function rect = makeOrientedRectangle(p1, p2, fraction, len, wid, offset)
    d = p2 - p1;
    d = d / norm(d);
    n = [-d(2), d(1)];
    c = p1 + fraction*(p2 - p1) + offset*n;
    rect = [c - 0.5*len*d - 0.5*wid*n;
            c + 0.5*len*d - 0.5*wid*n;
            c + 0.5*len*d + 0.5*wid*n;
            c - 0.5*len*d + 0.5*wid*n];
end
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
function plotPointWithLabel(p, labelText)
    plot(p(1), p(2), 'ko', 'MarkerFaceColor','w', 'MarkerSize',7, 'LineWidth',1.5);
    text(p(1)+10, p(2)+10, labelText, 'FontWeight','bold');
end
function drawCircle(center, radius, labelText)
    t = linspace(0, 2*pi, 100);
    x = center(1) + radius*cos(t);
    y = center(2) + radius*sin(t);
    plot(x, y, 'b--', 'LineWidth', 1.5);
    text(center(1)+radius+5, center(2), labelText, 'FontSize',8);
end
function exportFigure(figHandle, filename)
    try
        exportgraphics(figHandle, filename, 'Resolution', 300);
    catch
        saveas(figHandle, filename);
    end
end
