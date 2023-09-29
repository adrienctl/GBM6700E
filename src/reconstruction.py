import calibration
import numpy as np


def reconstruction_one_point(param_P0, param_P20, coord_P0,coord_P20):
    param_P0 = np.array(param_P0).flatten()
    param_P20 = np.array(param_P20).flatten()
    A = []
    B = []
    u_P0,v_P0 = coord_P0[0],coord_P0[1]
    u_P20,v_P20 = coord_P20[0],coord_P20[1]
    L1_P0,L2_P0, L3_P0, L4_P0, L5_P0, L6_P0, L7_P0, L8_P0, L9_P0, L10_P0, L11_P0 = param_P0[0],param_P0[1],param_P0[2],param_P0[3],param_P0[4],param_P0[5],param_P0[6],param_P0[7],param_P0[8],param_P0[9],param_P0[10]
    L1_P20,L2_P20, L3_P20, L4_P20, L5_P20, L6_P20, L7_P20, L8_P20, L9_P20, L10_P20, L11_P20 = param_P20[0],param_P20[1],param_P20[2],param_P20[3],param_P20[4],param_P20[5],param_P20[6],param_P20[7],param_P20[8],param_P20[9],param_P20[10]
    A.append([u_P0*L9_P0-L1_P0, u_P0*L10_P0-L2_P0, u_P0*L11_P0-L3_P0])
    A.append([v_P0*L9_P0-L5_P0, v_P0*L10_P0-L6_P0, v_P0*L11_P0-L7_P0])
    A.append([u_P20*L9_P20-L1_P20, u_P20*L10_P20-L2_P20, u_P20*L11_P20-L3_P20])
    A.append([v_P20*L9_P20-L5_P20, v_P20*L10_P20-L6_P20, v_P20*L11_P20-L7_P20])
    B.append([L4_P0-u_P0])
    B.append([L8_P0-v_P0])
    B.append([L4_P20-u_P20])
    B.append([L8_P20-v_P20])
    
    #Ici on a fait le rapport, on utilise pas certains param, je pense qu'on perd de l'info
    """
    A.append([coord_P0[0]*param_P0[4]-coord_P0[1]*param_P0[0],
            coord_P0[0]*param_P0[5]-coord_P0[1]*param_P0[1],
            coord_P0[0]*param_P0[6]-coord_P0[1]*param_P0[2]])
    B.append([param_P0[3]*coord_P0[1]-param_P0[7]*coord_P0[0]])
    A.append([coord_P20[0]*param_P20[4]-coord_P20[1]*param_P20[0],
            coord_P20[0]*param_P20[5]-coord_P20[1]*param_P20[1],
            coord_P20[0]*param_P20[6]-coord_P20[1]*param_P20[2]])
    B.append([param_P20[3]*coord_P20[1]-param_P20[7]*coord_P20[0]])
    
    
    A.append([param_P0[8]*(coord_P0[0]+coord_P0[1])-param_P0[1]-param_P0[4],
              param_P0[9]*(coord_P0[0]+coord_P0[1])-param_P0[2]-param_P0[5],
              param_P0[10]*(coord_P0[0]+coord_P0[1])-param_P0[3]-param_P0[6]])
    B.append([-coord_P0[0]-coord_P0[1]+param_P0[3]+param_P0[7]])
    A.append([param_P20[8]*(coord_P20[0]+coord_P20[1])-param_P20[1]-param_P20[4],
              param_P20[9]*(coord_P20[0]+coord_P20[1])-param_P20[2]-param_P20[5],
              param_P20[10]*(coord_P20[0]+coord_P20[1])-param_P20[3]-param_P20[6]])
    B.append([-coord_P20[0]-coord_P20[1]+param_P20[3]+param_P20[7]])

    A.append([param_P0[8]*(coord_P0[0]-coord_P0[1])-param_P0[1]+param_P0[4],
              param_P0[9]*(coord_P0[0]-coord_P0[1])-param_P0[2]+param_P0[5],
              param_P0[10]*(coord_P0[0]-coord_P0[1])-param_P0[3]+param_P0[6]])
    B.append([-coord_P0[0]+coord_P0[1]+param_P0[3]-param_P0[7]])
    A.append([param_P20[8]*(coord_P20[0]-coord_P20[1])-param_P20[1]+param_P20[4],
              param_P20[9]*(coord_P20[0]-coord_P20[1])-param_P20[2]+param_P20[5],
              param_P20[10]*(coord_P20[0]-coord_P20[1])-param_P20[3]+param_P20[6]])
    B.append([-coord_P20[0]+coord_P20[1]+param_P20[3]-param_P20[7]])
    """


    XYZ,_,_,_ = np.linalg.lstsq(A, B, rcond=None)
    return XYZ.flatten()

def reconstruction_all_beads(param_P0, param_P20,mat2D_P0,mat2D_P20):
    R = []
    for i in range(len(mat2D_P0)):
        for k in range(len(mat2D_P20)):
            if mat2D_P0[i][0] == mat2D_P20[k][0]: #tester si les noms sont les memes
                R.append(reconstruction_one_point(param_P0, param_P20,mat2D_P0[i][1][0],mat2D_P20[k][1][0]))
    return np.array(R)

def reconstruction_vertebrae(param_P0, param_P20,mat2D_P0,mat2D_P20):
    R = []
    for i in range(len(mat2D_P0)):
        for k in range(len(mat2D_P20)):
            if mat2D_P0[i][0][0] == mat2D_P20[k][0][0]: #tester si les noms sont les memes
                p_P0 = mat2D_P0[i][0][1][0]
                p_P20 = mat2D_P20[k][0][1][0]
                for j in range(6):
                    pts_P0 = p_P0[j][1][0]
                    pts_P20 = p_P20[j][1][0]
                    R.append(reconstruction_one_point(param_P0, param_P20,pts_P0,pts_P20))
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
