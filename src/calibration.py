import numpy as np
import random

def create_equations(mat2D,mat3D,nb_beads=np.inf):
    A = []
    B = []
    nb = 0

    random_ind = list(range(len(mat3D)))
    random.shuffle(random_ind)
    random_mat3D = [mat3D[i] for i in random_ind]
    

    #mat2D,mat3D = np.random.shuffle(mat2D),np.random.shuffle(mat3D)
    for line3D in random_mat3D:
        for line2D in mat2D:
            if line2D[0]==line3D[0] and nb<nb_beads:
                X,Y,Z = line3D[1][0][0],line3D[1][0][1],line3D[1][0][2]
                u,v = line2D[1][0][0],line2D[1][0][1]
                A.append([X,Y,Z,1,0,0,0,0,-u*X,-u*Y,-u*Z])
                B.append([u])
                A.append([0,0,0,0,X,Y,Z,1,-v*X,-v*Y,-v*Z])
                B.append([v])
                nb+=1
    return np.array(A),np.array(B)


    

def compute_camera_parameters(mat2D,mat3D,nb_beads=np.inf):
    """
    input: mat2D : image 2D des beads, mat3D : coordonnes 3D des beads
    output: M : matrice de paramètres de la caméra (contenant 11 paramètres intrasèques ou extrasèques)
    """
    A,B = create_equations(mat2D,mat3D,nb_beads)
    L,_,_,_ = np.linalg.lstsq(A, B, rcond=None)
    M = np.ones((3,4))
    for i in range(len(L)):
        M[int(i/4)][i%4] = L[i][0]
    return M


def select_config_beads(mat3D,config):
    """
    input: mat3D : image 3D des beads
    config : liste des 8 beads sélectionnés pour la calibration
    """
    mat3D_selected = []
    if config==None:
        return mat3D
    for line3D in mat3D:
        if line3D[0] in config:
            mat3D_selected.append(line3D)
    return mat3D_selected


def print_parameters(M):
    L = M.flatten()
    print("L =",L)
    L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11 = L[0],L[1],L[2],L[3],L[4],L[5],L[6],L[7],L[8],L[9],L[10]
    Axyz = np.linalg.pinv(np.array([[L1,L2,L3],[L5,L6,L7],[L9,L10,L11]]))
    Bxyz = np.array([-L4,-L8,-1])
    X0,Y0,Z0 = Axyz@Bxyz
    d = -1/(L9**2+L10**2+L11**2)**0.5
    u0 = d**2*(L1*L9+L2*L10+L3*L11)
    v0 = d**2*(L5*L9+L6*L10+L7*L11)
    cu = (d**2*((u0*L9-L1)**2+(u0*L10-L2)**2+(u0*L11-L3)**2))**0.5
    cv = (d**2*((v0*L9-L5)**2+(v0*L10-L6)**2+(v0*L11-L7)**2))**0.5
    R11 = (u0*L9-L1)*d/cu
    R12 = (u0*L10-L2)*d/cu
    R13 = (u0*L11-L3)*d/cu
    R21 = (v0*L9-L5)*d/cv
    R22 = (v0*L10-L6)*d/cv
    R23 = (v0*L11-L7)*d/cv
    R31 = d*L9
    R32 = d*L10
    R33 = d*L11
    R = np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]])

    print("____________________________________________________________________________")
    print("Paramètres intrinsèques de la caméra :")
    print("(u0,v0) = ({0:.2f},{1:.2f})".format(u0,v0))
    print("(cu,cv) = ({0:.2f},{1:.2f})".format(cu,cv))
    print("Paramètres extrasèques de la caméra :")
    print("(X0,Y0,Z0) = ("+str(X0)+","+str(Y0)+","+str(Z0)+")")
    print("Matrice de rotation R :", R)
    print("____________________________________________________________________________\n")

def add_2D_gaussian_noise(mat2D,sigma):
    for i in range(len(mat2D)):
        mat2D[i][1][0][0] += np.random.normal(0,sigma)
        mat2D[i][1][0][1] += np.random.normal(0,sigma)
    return mat2D

def add_3D_gaussian_noise(mat3D,sigma):
    for i in range(len(mat3D)):
        mat3D[i][1][0][0] += np.random.normal(0,sigma)
        mat3D[i][1][0][1] += np.random.normal(0,sigma)
        mat3D[i][1][0][2] += np.random.normal(0,sigma)
    return mat3D
