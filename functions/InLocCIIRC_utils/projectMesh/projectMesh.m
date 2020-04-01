function [RGBcut, XYZcut, depth] = projectMesh(meshPath, f, R, t, sensorSize, ortho, mag, projectMeshPyPath)
% TODO: support outputSize param. Then interpolation may be necessary for the XYZcut

inputPath = strcat(tempname, '.mat');
outputPath = strcat(tempname, '.mat');
save(inputPath, 'meshPath', 'f', 'R', 't', 'sensorSize', 'ortho', 'mag');

% call projectPointCloud.py
command = sprintf('PATH=/usr/local/bin:$PATH PYOPENGL_PLATFORM=osmesa python3 %s %s %s', projectMeshPyPath, inputPath, outputPath);
disp(command)
[status, cmdout] = system(command);
disp(cmdout)

% load results
load(outputPath, 'RGBcut', 'XYZcut', 'depth')

% delete temporary files
delete(inputPath);
delete(outputPath);

end
