import calibration
import reconstruction
import affichage
import lecture_ecriture
import errors
import time
import numpy as np
import matplotlib.pyplot as plt

config_1 = [["A_1_1","A_5_1","A_1_6","A_5_6","B_1_1","B_5_1","B_1_5","B_5_5"],"Config 1"]
config_2 = [["A_1_1","A_5_1","A_1_6","A_5_6","B_2_2","B_4_2","B_2_4","B_4_4"],"Config 2"]
config_3 = [["A_2_3","A_4_3","A_2_4","A_4_4","B_2_2","B_4_2","B_2_4","B_4_4"],"Config 3"]
config_4 = [["A_1_1","A_2_1","A_1_2","A_2_2","B_4_5","B_5_4","B_4_4","B_5_5"],"Config 4"]
config_5 = [["A_1_1","A_5_1","A_1_6","A_5_6","A_2_3","A_4_3","A_2_4","A_4_4"],"Config 5"]
config_6 = [["B_2_2","B_4_2","B_2_4","B_4_4","B_1_1","B_5_1","B_1_5","B_5_5"],"Config 6"]
config_0 = [None,"All beads"]

CONFIG = config_1
SUBPLOT = True # pour faire des plots individuels pour le rapport

def main():
    print("### Début de la calibration des cameras ###\n")
    start_time = time.time()
    nb_beads_calibration_max = 54
    Beads2D_calib_PA0, Beads2D_calib_PA20 = lecture_ecriture.load_calib_2D("data/Calib_Beads2D.mat")
    Beads3D_calib = lecture_ecriture.load_calib_3D("data/Calib_Beads3D.mat")
    Beads3D_calib_selected = calibration.select_config_beads(Beads3D_calib,CONFIG[0])
    param_camera_PA0 = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib_selected,nb_beads_calibration_max)
    param_camera_PA20 = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib_selected,nb_beads_calibration_max)
    param_camera_PA0_groundtruth = calibration.compute_camera_parameters(Beads2D_calib_PA0,Beads3D_calib)
    param_camera_PA20_groundtruth = calibration.compute_camera_parameters(Beads2D_calib_PA20,Beads3D_calib)
    print("\nCAMERA PA0 :")
    calibration.print_parameters(param_camera_PA0)
    print("\nCAMERA PA20 :")
    calibration.print_parameters(param_camera_PA20)
    fin_calib_time = time.time()
    print("Fin de la calibration des cameras, éxécutée en {0:.2f} secondes.\n".format(fin_calib_time-start_time))

    print("### Début de la reconstruction 3D ###\n")
    Beads2D_vert_PA0, Beads2D_vert_PA20 = lecture_ecriture.load_vert_2D("data/Vertebrae2D.mat")
    vert_3D = reconstruction.reconstruction_vertebrae(param_camera_PA0,param_camera_PA20,Beads2D_vert_PA0,Beads2D_vert_PA20)
    vert_3D_groundtruth = reconstruction.reconstruction_vertebrae(param_camera_PA0_groundtruth,param_camera_PA20_groundtruth,Beads2D_vert_PA0,Beads2D_vert_PA20)
    # test to reconstruct the beads
    calib_3D_test = reconstruction.reconstruction_all_beads(param_camera_PA0,param_camera_PA20,Beads2D_calib_PA0,Beads2D_calib_PA20)
    fin_reconstr_time = time.time()
    print("Fin de la reconstruction 3D, éxécutée en {0:.2f} secondes.\n".format(fin_reconstr_time-fin_calib_time))

    print("### Début de l'affichage' ###\n")
    affichage.plot_3D_points(vert_3D, vert_3D_groundtruth,Beads3D_calib,"Vertebrae reconstruction with "+CONFIG[1],CONFIG[0],SUBPLOT)
    affichage.create_RMS_curve(SUBPLOT)
    affichage.create_2D_noise_curve(SUBPLOT)
    affichage.create_3D_noise_curve(SUBPLOT)
    affichage.plot_errors_bary(Beads3D_calib_selected, vert_3D, vert_3D_groundtruth,SUBPLOT)
    if SUBPLOT :
        manager = plt.get_current_fig_manager() # Pour mettre la fenêtre en plein écran
        manager.full_screen_toggle()
        plt.show()
    print("Fin de l'affichage utilisé pendant {0:.2f} secondes.\n".format(time.time()-fin_reconstr_time))



if __name__ == "__main__":
    main()
    


