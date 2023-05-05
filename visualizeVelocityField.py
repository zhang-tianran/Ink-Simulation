import numpy as np
import matplotlib.pyplot as plt

# TODO: Change this to read in grid dimensions instead of hard-coding
WATERGRID_X = 2
WATERGRID_Y = 2
WATERGRID_Z = 2


# Returns two numpy arrays:
#   1) The starting positions of your velocity vectors
#   2) The directions of the velocity vectors
def load_water_grid():
    positions = np.array([
        (0, 0, 0), \
        (1, 0, 0), \
        (0, 1, 0), \
        (0, 0, 1), \
        (1, 1, 0), \
        (0, 1, 1), \
        (1, 0, 1), \
        (1, 1, 1), \
        ])
    directions = np.array([
        (1, 1, 1), \
        (1, 1, 1), \
        (1, 1, 1), \
        (1, 1, 1), \
        (1, 1, 1), \
        (1, 1, 1), \
        (1, 1, 1), \
        (1, 1, 1), \
        ])
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