import numpy as np
import matplotlib.pyplot as plt
from projectMesh import projectMesh

meshPath = '/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/models/B-670/mesh_rotated.obj'
f = 600
R = np.eye(3)
t = np.transpose(np.array([0.0, 1.0, 0.0]))
sensorSize = np.array([1600, 1200])

RGBcut, XYZcut, depth = projectMesh(meshPath, f, R, t, sensorSize, False, -1)

plt.figure()
plt.imshow(RGBcut)

plt.figure()
plt.imshow(depth)

plt.figure()
plt.title('Orthographic projection')
RGBcut, XYZcut, depth = projectMesh(meshPath, f, R, t, sensorSize, True, 3.0)
plt.imshow(RGBcut)

plt.show()