import matplotlib.pyplot as plt
import numpy as np
import errors

def plot_3D_points(mat3D_vert, mat3D_vert_groundtruth,mat3D_beads,name,config):
    fig = plt.figure(name)
    ax = fig.add_subplot(111, projection='3d')
    X,Y,Z = mat3D_vert.T[0],mat3D_vert.T[1],mat3D_vert.T[2]
    ax.scatter(X,Y,Z,label = "Vertebrae reconstruction",c='b',marker='o')
    X,Y,Z = mat3D_vert_groundtruth.T[0],mat3D_vert_groundtruth.T[1],mat3D_vert_groundtruth.T[2]
    ax.scatter(X,Y,Z,label = "Groundtruth vertebrae reconstruction",c='b',marker='+',alpha=0.25)
    for line in mat3D_beads:
        if config!=None and line[0] in config:
            ax.scatter(line[1][0][0],line[1][0][1],line[1][0][2],c='g',marker='o')
        else :
            ax.scatter(line[1][0][0],line[1][0][1],line[1][0][2],c='r',marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    max_range = np.array([mat3D_vert[:,0].max()-mat3D_vert[:,0].min(), mat3D_vert[:,1].max()-mat3D_vert[:,1].min(), mat3D_vert[:,2].max()-mat3D_vert[:,2].min()]).max() /1.5
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    #manager = plt.get_current_fig_manager() # Pour mettre la fenêtre en plein écran
    #manager.full_screen_toggle()
    plt.legend()
    #plt.show()

def plot_errors_bary(beads3D_selected, vert_3D, vert_3D_groundtruth):
    plt.subplot(1,3,3)
    dist, errors_list = errors.dist_bary_gt(beads3D_selected, vert_3D, vert_3D_groundtruth)
    plt.plot(dist, errors_list, 'ro')
    plt.xlabel("Distance to barycenter")
    plt.ylabel("Error")
    plt.title("Error depending on distance to barycenter")
    #plt.show()