import calibration
import numpy as np


def reconstruction_one_point(param_P0, param_P20, coord_P0,coord_P20):
    param_P0 = np.array(param_P0).flatten()
    param_P20 = np.array(param_P20).flatten()
    A = []
    B = []
    A.append([coord_P0[0]*param_P0[4]-coord_P0[1]*param_P0[0],
            coord_P0[0]*param_P0[5]-coord_P0[1]*param_P0[1],
            coord_P0[0]*param_P0[6]-coord_P0[1]*param_P0[2]])
    B.append([param_P0[3]*coord_P0[1]-param_P0[2]*coord_P0[0]])
    A.append([coord_P20[0]*param_P20[4]-coord_P20[1]*param_P20[0],
            coord_P20[0]*param_P20[5]-coord_P20[1]*param_P20[1],
            coord_P20[0]*param_P20[6]-coord_P20[1]*param_P20[2]])
    B.append([param_P20[3]*coord_P20[1]-param_P20[2]*coord_P20[0]])
    XYZ,_,_,_ = np.linalg.lstsq(A, B, rcond=None)
    return np.array(XYZ)

def reconstruction_all_beads(param_P0, param_P20,mat2D_P0,mat2D_P20):
    R = []
    for i in range(len(mat2D_P0)):
        R.append(reconstruction_one_point(param_P0, param_P20,mat2D_P0[i][1][0],mat2D_P20[i][1][0]))
    return np.array(R)

def reconstruction_vertebrae(param_P0, param_P20,mat2D_P0,mat2D_P20):
    R = []
    for i in range(len(mat2D_P0)):
        for j in range(6):
            p_P0 = mat2D_P0[i][0][1][0][j]
            p_P20 = mat2D_P20[i][0][1][0][j]
            R.append(reconstruction_one_point(param_P0, param_P20,p_P0[1][0],p_P20[1][0]))
    return np.array(R)




def reconstruction_all_beads_old(mat2D, M):
    R = []
    for coord2D in mat2D:
        R.append(reconstruction_one_point([coord2D[1][0][0],coord2D[1][0][1],1],M))
    return np.array(R)

def reconstruction_vertebrae_old(mat2D,M):
    L = []
    for vertebre in mat2D:
        for point in vertebre[0][1][0]:
            L.append(reconstruction_one_point([point[1][0][0],point[1][0][1],1],M))
    return np.array(L)
