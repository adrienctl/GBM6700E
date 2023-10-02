import matplotlib.pyplot as plt
import numpy as np
import errors

def plot_3D_points(mat3D_vert, mat3D_vert_groundtruth,mat3D_beads,name,config,SUBPLOT):
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
    ax.set_title(name)
    ax.set_xlim(-600, 0)
    ax.set_ylim(-300, 300)
    ax.set_zlim(-300, 300)
    ax.grid(False)
    #manager = plt.get_current_fig_manager() # Pour mettre la fenêtre en plein écran
    #manager.full_screen_toggle()
    plt.legend()
    if not SUBPLOT : plt.show()
    #plt.show()

def plot_errors_bary(beads3D_selected, vert_3D, vert_3D_groundtruth,SUBPLOT):
    if SUBPLOT : plt.subplot(3,4,4)
    dist, errors_list_3D, errors_list_X, errors_list_Y, errors_list_Z = errors.dist_bary_gt(beads3D_selected, vert_3D, vert_3D_groundtruth)
    plt.scatter(dist, errors_list_3D, label="3D error",s=3)
    plt.scatter(dist, errors_list_X, label="X error",s=3)
    plt.scatter(dist, errors_list_Y, label="Y error",s=3)
    plt.scatter(dist, errors_list_Z, label="Z error",s=3)
    plt.xlabel("Distance to barycenter")
    plt.ylabel("Error")
    plt.title("Error depending on distance to barycenter")
    plt.legend()
    if not SUBPLOT : plt.show()
    #plt.show()