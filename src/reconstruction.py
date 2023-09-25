import calibration
import numpy as np


def reconstruction_one_point(coord2D,M):
    inv_M = (np.linalg.inv(M.T.dot(M))).dot(M.T)
    return np.dot(inv_M,coord2D)

def reconstruction_all_beads(mat2D, M):
    L = []
    for coord2D in mat2D:
        L.append(reconstruction_one_point([coord2D[1][0][0],coord2D[1][0][1],1],M))
    return np.array(L)

def reconstruction_vertebrae(mat2D,M):
    L = []
    for vertebre in mat2D:
        for point in vertebre[0][1][0]:
            L.append(reconstruction_one_point([point[1][0][0],point[1][0][1],1],M))
    return np.array(L)
