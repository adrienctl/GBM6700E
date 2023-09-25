import matplotlib.pyplot as plt

def plot_3D_points(mat3D,name):
    fig = plt.figure(name)
    ax = fig.add_subplot(111, projection='3d')
    for line in mat3D:
        ax.scatter(line[0],line[1],line[2],c='b',marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    max_range = max(max(line) for line in mat3D)
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])
    plt.show()