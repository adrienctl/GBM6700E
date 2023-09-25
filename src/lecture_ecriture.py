import scipy.io as sio

def load_calib_2D(file_name):
    mat_contents = sio.loadmat(file_name)
    Beads2D_calib_PA0 = mat_contents['Beads2D_PA0'][0]
    Beads2D_calib_PA20 = mat_contents['Beads2D_PA20'][0]
    return Beads2D_calib_PA0, Beads2D_calib_PA20

def load_calib_3D(file_name):
    mat_contents = sio.loadmat(file_name)
    Beads3D_calib = mat_contents['Calib_Beads3D'][0]
    return Beads3D_calib

def load_vert_2D(file_name):
    mat_contents = sio.loadmat(file_name)
    Beads2D_vert_PA0 = mat_contents['Vertebrae_PA0']
    Beads2D_vert_PA20 = mat_contents['Vertebrae_PA20']
    return Beads2D_vert_PA0, Beads2D_vert_PA20


