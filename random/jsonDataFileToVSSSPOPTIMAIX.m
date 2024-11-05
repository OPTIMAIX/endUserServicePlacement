clear;

% Load JSON data from a file
filename = 'input_1_3_3.json';
json_data = fileread(filename);

placement = NFPlanningSolverStructs(json_data);

