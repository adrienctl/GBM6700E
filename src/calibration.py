import numpy as np

def create_equations(mat2D,mat3D):
    A = []
    B = []
    for line3D in mat3D:
        for line2D in mat2D:
            if line2D[0]==line3D[0]:
                A.append([line3D[1][0][0],
                          line3D[1][0][1],
                          line3D[1][0][2],
                          1,
                          0,
                          0,
                          0,
                          0,
                          -line2D[1][0][0]*line3D[1][0][0],
                          -line2D[1][0][0]*line3D[1][0][1],
                          -line2D[1][0][0]*line3D[1][0][2]])
                B.append([line2D[1][0][0]])
                A.append([0,
                          0,
                          0,
                          0,
                          line3D[1][0][0],
                          line3D[1][0][1],
                          line3D[1][0][2],
                          1,
                          -line2D[1][0][1]*line3D[1][0][0],
                          -line2D[1][0][1]*line3D[1][0][1],
                          -line2D[1][0][1]*line3D[1][0][2]])
                B.append([line2D[1][0][1]])
    return np.array(A),np.array(B)
    

def common_beads(mat2D_PA0, mat2d_PA20):
    common_beads_list=[]
    for line0 in mat2D_PA0:
        for line20 in mat2d_PA20:
            if line0[0]==line20[0]:
                common_beads_list.append(line0[0])
    return common_beads_list

def compute_camera_parameters(mat2D,mat3D):
    """
    input: mat2D : image 2D des beads, mat3D : coordonnes 3D des beads
    output: M : matrice de paramètres de la caméra (contenant 11 paramètres intrasèques ou extrasèques)
    """
    A,B = create_equations(mat2D,mat3D)
    L,_,_,_ = np.linalg.lstsq(A, B, rcond=None)
    M = np.ones((3,4))
    for i in range(len(L)):
        M[int(i/4)][i%4] = L[i][0]
    return M

def print_parameters(M):
    L = M.flatten()
    d = -1/(L[8]**2+L[9]**2+L[10]**2)**0.5
    u0 = d**2*(L[0]*L[8]+L[1]*L[9]+L[2]*L[10])
    v0 = d**2*(L[4]*L[8]+L[5]*L[9]+L[6]*L[10])
    cu = (d**2*((u0*L[8]-L[0])**2+(u0*L[9]-L[1])**2+(u0*L[10]-L[2])**2))**0.5
    cv = (d**2*((v0*L[8]-L[4])**2+(v0*L[9]-L[5])**2+(v0*L[10]-L[6])**2))**0.5
    Axyz = np.array([[L[0],L[1],L[2]],[L[4],L[5],L[6]],[L[8],L[9],L[10]]])
    Bxyz = np.array([[L[3]],[L[7]],[1]])
    X0,Y0,Z0 = np.dot(np.linalg.inv(Axyz),Bxyz)
    R11 = (u0*L[8]-L[0])*d/cu
    R12 = (u0*L[9]-L[1])*d/cu
    R13 = (u0*L[10]-L[2])*d/cu
    R21 = (v0*L[8]-L[4])*d/cv
    R22 = (v0*L[9]-L[5])*d/cv
    R23 = (v0*L[10]-L[6])*d/cv
    R31 = d*L[8]
    R32 = d*L[9]
    R33 = d*L[10]
    R = np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]])


    print("____________________________________________________________________________")
    print("Paramètres intrinsèques de la caméra :")
    print("(u0,v0) = ({0:.2f},{1:.2f})".format(u0,v0))
    print("(cu,cv) = ({0:.2f},{1:.2f})".format(cu,cv))
    print("Paramètres extrasèques de la caméra :")
    print("(X0,Y0,Z0) = ("+str(X0[0])+","+str(Y0[0])+","+str(Z0[0])+")")
    print("Matrice de rotation R :", R)
    print("____________________________________________________________________________\n")