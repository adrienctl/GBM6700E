import calibration
import reconstruction
import affichage
import lecture_ecriture
import time

def main():
    print("### Début de la calibration des cameras ###\n")
    start_time = time.time()

    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib)

    print("\nParamètres de la caméra PA0 :")
    calibration.print_parameters(calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib))
    print("\nParamètres de la caméra PA20 :")
    calibration.print_parameters(calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib))
    fin_calib_time = time.time()
    print("Fin de la calibration des cameras, éxécutée en {0:.2f} secondes.\n".format(fin_calib_time-start_time))
    print("### Début de la reconstruction 3D ###\n")

    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")
    vert_3D = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)

    fin_reconstr_time = time.time()
    print("Fin de la reconstruction 3D, éxécutée en {0:.2f} secondes.\n".format(fin_reconstr_time-fin_calib_time))
    affichage.plot_3D_points(vert_3D)

    print("Utilisation de l'affichage pendant {0:.2f} secondes.\n".format(time.time()-fin_reconstr_time))


if __name__ == "__main__":
    main()
    