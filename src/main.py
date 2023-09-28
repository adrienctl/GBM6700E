import calibration
import reconstruction
import affichage
import lecture_ecriture
import errors
import time
import numpy as np
import matplotlib.pyplot as plt

def create_RMS_curve():
    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")

    # calcul du groundtruth
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib,50)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib,50)
    vert_3D_groundtruth = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)

    print(vert_3D_groundtruth.shape)

    RMS_X = []
    RMS_Y = []
    RMS_Z = []
    RMS_3D = []

    abs = np.arange(1,50,1)

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
    plt.legend()
    plt.show()

def main():
    print("### Début de la calibration des cameras ###\n")
    start_time = time.time()

    nb_beads_calibration_max = 17
    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib,nb_beads_calibration_max)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib,nb_beads_calibration_max)

    print("\nParamètres de la caméra PA0 :")
    calibration.print_parameters(param_camera_PA0)
    print("\nParamètres de la caméra PA20 :")
    calibration.print_parameters(param_camera_PA20)
    fin_calib_time = time.time()
    print("Fin de la calibration des cameras, éxécutée en {0:.2f} secondes.\n".format(fin_calib_time-start_time))
    print("### Début de la reconstruction 3D ###\n")

    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")
    vert_3D = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)
    calib_3D_test = reconstruction.reconstruction_all_beads(param_camera_PA0,param_camera_PA20,Beads2D_calib_PA0,Beads2D_calib_PA20)
    fin_reconstr_time = time.time()
    print("Fin de la reconstruction 3D, éxécutée en {0:.2f} secondes.\n".format(fin_reconstr_time-fin_calib_time))
    affichage.plot_3D_points(vert_3D,"Vertebrae reconstruction")
    affichage.plot_3D_points(calib_3D_test,"Beads reconstruction")
    end_display_time = time.time()

    print("Utilisation de l'affichage pendant {0:.2f} secondes.\n".format(end_display_time-fin_reconstr_time))

    create_RMS_curve()

    print("Fin du programme.")


if __name__ == "__main__":
    main()
    


