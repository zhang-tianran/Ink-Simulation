import numpy as np
import matplotlib.pyplot as plt
from plyfile import PlyData, PlyElement

# Returns two numpy arrays:
#   1) The starting positions of your velocity vectors
#   2) The directions of the velocity vectors
def load_water_grid():
    plydata = PlyData.read('./output/velocityField.ply')

    x = np.array(plydata['vertex']['x'])
    y = np.array(plydata['vertex']['y'])
    z = np.array(plydata['vertex']['z'])
    positions = np.dstack((x, y, z))[0]

    u = np.array(plydata['vertex']['u'])
    v = np.array(plydata['vertex']['v'])
    w = np.array(plydata['vertex']['w'])
    directions = np.dstack((u, v, w))[0]
    return positions, directions


def render_velocity_field():
    positions, directions = load_water_grid()
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(positions[:, 0], positions[:, 1], positions[:, 2], \
              directions[:, 0], directions[:, 1], directions[:, 2], \
              length=0.2, normalize=False)
    ax.set_title("Water Grid Velocity Field")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()


if __name__ == "__main__":
    render_velocity_field()