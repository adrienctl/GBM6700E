import matplotlib.pyplot as plt
import numpy as np
import errors
import calibration
import reconstruction
import lecture_ecriture

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

def create_RMS_curve(SUBPLOT):
    if SUBPLOT : plt.subplot(3,4,1)
    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")

    # calcul du groundtruth
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib,50)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib,50)
    vert_3D_groundtruth = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)

    RMS_X = []
    RMS_Y = []
    RMS_Z = []
    RMS_3D = []

    abs = np.arange(6,51,1)

    for nb_beads in abs:
        param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib,nb_beads)
        param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib,nb_beads)
        vert_3D = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)
        RMS_X.append(errors.RMS(vert_3D,vert_3D_groundtruth,0))
        RMS_Y.append(errors.RMS(vert_3D,vert_3D_groundtruth,1))
        RMS_Z.append(errors.RMS(vert_3D,vert_3D_groundtruth,2))
        RMS_3D.append(errors.RMS(vert_3D,vert_3D_groundtruth))
    
    plt.plot(abs,RMS_X,label="RMS X",color="red", linestyle=':')
    plt.plot(abs,RMS_Y,label="RMS Y",color="green", linestyle=':')
    plt.plot(abs,RMS_Z,label="RMS Z",color="blue", linestyle=':')
    plt.plot(abs,RMS_3D,label="RMS 3D",color="black", linestyle='solid')
    plt.xlabel("Nombre de beads utilisés pour la calibration")
    plt.ylabel("RMS")
    plt.title("RMS en fonction du nombre de beads utilisés")
    plt.legend()
    if not SUBPLOT : plt.show()



def create_2D_noise_curve(SUBPLOT):
    sigma_max = 0.5

    if SUBPLOT : plt.subplot(3,4,5)
    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")
    # calcul du groundtruth
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib)
    vert_3D_groundtruth = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)

    sigma_list = []
    RMS_X = []
    RMS_Y = []
    RMS_Z = []
    RMS_3D = []
    for sigma in np.linspace(0,sigma_max,50):
        sigma_list.append(sigma)
        Beads2D_calib_PA0_noise, Beads2D_calib_PA20_noise = calibration.add_2D_gaussian_noise(Beads2D_calib_PA0,sigma), calibration.add_2D_gaussian_noise(Beads2D_calib_PA20,sigma)
        param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0_noise,Beads3D_calib)
        param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20_noise,Beads3D_calib)
        vert_3D = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)
        RMS_X.append(errors.RMS(vert_3D,vert_3D_groundtruth,0))
        RMS_Y.append(errors.RMS(vert_3D,vert_3D_groundtruth,1))
        RMS_Z.append(errors.RMS(vert_3D,vert_3D_groundtruth,2))
        RMS_3D.append(errors.RMS(vert_3D,vert_3D_groundtruth))
    
    plt.plot(sigma_list,RMS_X,label="RMS X",color="red", linestyle=':')
    plt.plot(sigma_list,RMS_Y,label="RMS Y",color="green", linestyle=':')
    plt.plot(sigma_list,RMS_Z,label="RMS Z",color="blue", linestyle=':')
    plt.plot(sigma_list,RMS_3D,label="RMS 3D",color="black", linestyle='solid')
    plt.xlabel("Sigma")
    plt.ylabel("Error")
    plt.title("Effect of noise on 2D points on reconstruction")
    plt.legend()
    if not SUBPLOT : plt.show()
    # with 8 beads 

    if SUBPLOT : plt.subplot(3,4,9)
    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")
    # calcul du groundtruth
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib)
    vert_3D_groundtruth = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)

    sigma_list = []
    RMS_X = []
    RMS_Y = []
    RMS_Z = []
    RMS_3D = []
    for sigma in np.linspace(0,sigma_max,50):
        sigma_list.append(sigma)
        Beads2D_calib_PA0_noise, Beads2D_calib_PA20_noise = calibration.add_2D_gaussian_noise(Beads2D_calib_PA0,sigma), calibration.add_2D_gaussian_noise(Beads2D_calib_PA20,sigma)
        param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0_noise,Beads3D_calib,8)
        param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20_noise,Beads3D_calib,8)
        vert_3D = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)
        RMS_X.append(errors.RMS(vert_3D,vert_3D_groundtruth,0))
        RMS_Y.append(errors.RMS(vert_3D,vert_3D_groundtruth,1))
        RMS_Z.append(errors.RMS(vert_3D,vert_3D_groundtruth,2))
        RMS_3D.append(errors.RMS(vert_3D,vert_3D_groundtruth))
    
    plt.plot(sigma_list,RMS_X,label="RMS X",color="red", linestyle=':')
    plt.plot(sigma_list,RMS_Y,label="RMS Y",color="green", linestyle=':')
    plt.plot(sigma_list,RMS_Z,label="RMS Z",color="blue", linestyle=':')
    plt.plot(sigma_list,RMS_3D,label="RMS 3D",color="black", linestyle='solid')
    plt.xlabel("Sigma")
    plt.ylabel("Error")
    plt.title("Effect of noise on 2D with 8 beads")
    plt.legend()
    if not SUBPLOT : plt.show()


def create_3D_noise_curve(SUBPLOT):
    sigma_max = 0.5

    if SUBPLOT : plt.subplot(3,4,8)
    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")
    # calcul du groundtruth
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib)
    vert_3D_groundtruth = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)

    sigma_list = []
    RMS_X = []
    RMS_Y = []
    RMS_Z = []
    RMS_3D = []
    for sigma in np.linspace(0,sigma_max,50):
        sigma_list.append(sigma)
        Beads3D_calib_noise = calibration.add_3D_gaussian_noise(Beads3D_calib,sigma)
        param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib_noise)
        param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib_noise)
        vert_3D = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)
        RMS_X.append(errors.RMS(vert_3D,vert_3D_groundtruth,0))
        RMS_Y.append(errors.RMS(vert_3D,vert_3D_groundtruth,1))
        RMS_Z.append(errors.RMS(vert_3D,vert_3D_groundtruth,2))
        RMS_3D.append(errors.RMS(vert_3D,vert_3D_groundtruth))
    
    plt.plot(sigma_list,RMS_X,label="RMS X",color="red", linestyle=':')
    plt.plot(sigma_list,RMS_Y,label="RMS Y",color="green", linestyle=':')
    plt.plot(sigma_list,RMS_Z,label="RMS Z",color="blue", linestyle=':')
    plt.plot(sigma_list,RMS_3D,label="RMS 3D",color="black", linestyle='solid')
    plt.xlabel("Sigma")
    plt.ylabel("Error")
    plt.title("Effect of noise on 3D points on reconstruction")
    plt.legend()
    if not SUBPLOT : plt.show()

    #with 8 beads

    if SUBPLOT : plt.subplot(3,4,12)
    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")
    # calcul du groundtruth
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib)
    vert_3D_groundtruth = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)

    sigma_list = []
    RMS_X = []
    RMS_Y = []
    RMS_Z = []
    RMS_3D = []
    for sigma in np.linspace(0,sigma_max,50):
        sigma_list.append(sigma)
        Beads3D_calib_noise = calibration.add_3D_gaussian_noise(Beads3D_calib,sigma)
        param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib_noise,8)
        param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib_noise,8)
        vert_3D = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)
        RMS_X.append(errors.RMS(vert_3D,vert_3D_groundtruth,0))
        RMS_Y.append(errors.RMS(vert_3D,vert_3D_groundtruth,1))
        RMS_Z.append(errors.RMS(vert_3D,vert_3D_groundtruth,2))
        RMS_3D.append(errors.RMS(vert_3D,vert_3D_groundtruth))
    
    plt.plot(sigma_list,RMS_X,label="RMS X",color="red", linestyle=':')
    plt.plot(sigma_list,RMS_Y,label="RMS Y",color="green", linestyle=':')
    plt.plot(sigma_list,RMS_Z,label="RMS Z",color="blue", linestyle=':')
    plt.plot(sigma_list,RMS_3D,label="RMS 3D",color="black", linestyle='solid')
    plt.xlabel("Sigma")
    plt.ylabel("Error")
    plt.title("Effect of noise on 3D with 8 beads")
    plt.legend()
    if not SUBPLOT : plt.show()