import numpy as np

def RMS(test, groundtruth, nb_dim=3):
    """
    Compute the Root Mean Square error between two 1D arrays.
    nb_dim = 0 pour X, 1 pour Y, 2 pour Z
    """
    assert len(test) == len(groundtruth), "Les deux tableaux doivent avoir la même taille"
    if (nb_dim == 0):
        Xtest = test.T[0]
        XGT = groundtruth.T[0]
        return np.sqrt(np.sum((Xtest-XGT)**2)/len(Xtest))
    elif (nb_dim == 1):
        Ytest = test.T[1]
        YGT = groundtruth.T[1]
        return np.sqrt(np.sum((Ytest-YGT)**2)/len(Ytest))
    elif (nb_dim == 2):
        Ztest = test.T[2]
        ZGT = groundtruth.T[2]
        return np.sqrt(np.sum((Ztest-ZGT)**2)/len(Ztest))
    else:
        return np.sqrt(np.sum((test-groundtruth)**2)/len(test))
    



def dist_bary_gt(beads3D_selected, vert_3D, vert_3D_groundtruth):
    dist = []
    errors_list_3D = []
    errors_list_X = []
    errors_list_Y = []
    errors_list_Z = []
    bead_list = []
    for i in range(len(beads3D_selected)):
        bead_list.append(beads3D_selected[i][1][0])
    bary_beads = np.mean(bead_list, axis=0)
    for i in range(len(vert_3D)):
        dist.append(np.linalg.norm(vert_3D[i] - bary_beads))
        errors_list_3D.append(np.linalg.norm(vert_3D[i] - vert_3D_groundtruth[i]))
        errors_list_X.append(np.linalg.norm(vert_3D[i][0] - vert_3D_groundtruth[i][0]))
        errors_list_Y.append(np.linalg.norm(vert_3D[i][1] - vert_3D_groundtruth[i][1]))
        errors_list_Z.append(np.linalg.norm(vert_3D[i][2] - vert_3D_groundtruth[i][2]))
    return dist, errors_list_3D, errors_list_X, errors_list_Y, errors_list_Z


    
    
    