clc;
close all;
cd(fileparts(mfilename('fullpath')));
disp('REDESIGN AND VALIDATION PIPELINE STARTED....');
localTopOpt = fullfile(pwd, 'ebike_topology_optimization_outputs');
sourceTopOpt = fullfile(pwd, '..', 'topology_optimization', 'ebike_topology_optimization_outputs');
if ~exist(fullfile(localTopOpt, 'ebike_topology_optimization_step6.mat'), 'file')
    if exist(fullfile(sourceTopOpt, 'ebike_topology_optimization_step6.mat'), 'file')
        if exist(localTopOpt, 'dir'); rmdir(localTopOpt, 's'); end
        copyfile(sourceTopOpt, localTopOpt);
        disp('Copied topology input from ../topology_optimization/ebike_topology_optimization_outputs.');
    else
        error(['Cannot find topology optimization input. Run ../topology_optimization/main.m first, ', ...
               'or provide ebike_topology_optimization_outputs/ebike_topology_optimization_step6.mat.']);
    end
end
disp('Step 1: Interpret topology density into the final clean redesign');
run('interpret_topology_redesign.m');
disp('Step 2: Validate the final clean redesign using FEM');
run('validate_redesign_fem.m');
disp('REDESIGN AND VALIDATION PIPELINE COMPLETED. Final design = First Clean Redesign.');
