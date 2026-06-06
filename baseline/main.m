clc; 
close all;

cd(fileparts(mfilename('fullpath')));

disp(' BASELINE PIPELINE STARTED .....');
disp('Step 1: Build reference skeleton and topology design domain');
run('build_geometry_model.m');

disp('Step  2: Generate structured Q4 FE mesh');
run('generate_q4_mesh.m');

disp('Step 3: Run baseline FEM');
run('run_baseline_fem.m');

disp(' BASELINE PIPELINE COMPLETED ');
