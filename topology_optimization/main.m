clc; 
close all;

cd(fileparts(mfilename('fullpath')));
disp('TOPOLOGY OPTIMIZATION PIPELINE STARTED ....');
localBaseline = fullfile(pwd, 'ebike_baseline_fem_outputs');
sourceBaseline = fullfile(pwd, '..', 'baseline', 'ebike_baseline_fem_outputs');
if ~exist(fullfile(localBaseline, 'ebike_baseline_fem_step4_5.mat'), 'file')
    if exist(fullfile(sourceBaseline, 'ebike_baseline_fem_step4_5.mat'), 'file')
        if exist(localBaseline, 'dir'); rmdir(localBaseline, 's'); end
        copyfile(sourceBaseline, localBaseline);
        disp('Copied baseline input from ../baseline/ebike_baseline_fem_outputs.');
    else
        error(['Cannot find baseline input. Run ../baseline/main.m first, ', ...
               'or provide ebike_baseline_fem_outputs/ebike_baseline_fem_step4_5.mat.']);
    end
end
disp('Stage 4: Run SIMP topology optimization');
run('run_topology_optimization.m');
disp(' TOPOLOGY OPTIMIZATION PIPELINE COMPLETED');
