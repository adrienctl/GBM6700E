import matplotlib.pyplot as plt
import numpy as np

def plot_3D_points(mat3D_vert,mat3D_beads,name,selected_beads = []):
    fig = plt.figure(name)
    ax = fig.add_subplot(111, projection='3d')
    for line in mat3D_vert:
        #print(line[0])
        ax.scatter(line[0],line[1],line[2],c='b',marker='o')
    for line in mat3D_beads:
        ax.scatter(line[0],line[1],line[2],c='r',marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    max_range = np.array([mat3D_vert[:,0].max()-mat3D_vert[:,0].min(), mat3D_vert[:,1].max()-mat3D_vert[:,1].min(), mat3D_vert[:,2].max()-mat3D_vert[:,2].min()]).max() /1.5
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    plt.show()