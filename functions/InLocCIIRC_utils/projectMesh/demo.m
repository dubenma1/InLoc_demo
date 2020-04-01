meshPath = '/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/models/B-670/mesh_rotated.obj';
f = 600;
R = eye(3);
t = [0.0; 1.0; 0.0];
sensorSize = [1600 1200];
projectMeshPyPath = './projectMesh.py';

[RGBcut, XYZcut, depth] = projectMesh(meshPath, f, R, t, sensorSize, false, -1, projectMeshPyPath);

figure(1);
imshow(RGBcut);

figure(2);
imagesc(depth);