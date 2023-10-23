# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 17:34:53 2023

@author: victo
"""
#%% Imports 
import matplotlib.pyplot as plt
import numpy as np 
import csv
import random
import copy
from tqdm import tqdm
#%% Variables 
root = r'C:\CADOURS_VICTOR\Etudes\3A\3D reconstruction\data_beads_vertebres'

#%% fonctions pour lecture des fichiers
def lecture_fichier_2D(chemin):
    dictionnaire = {}
    with open (chemin, newline ='') as fichier:
        csvreader = csv.reader(fichier,delimiter=';')
        for row in csvreader: 
            cle = row[0].removeprefix("'").removesuffix("'")
            valeur_int = row[1].removeprefix("[").removesuffix("]")
            valeur_liste = valeur_int.split(',')
            valeur = [float(i) for i in valeur_liste]
            dictionnaire[cle] = valeur
    return dictionnaire

def lecture_fichier_3D(chemin):
    dictionnaire = {}
    liste_vertebres =[]
    liste_6_points = ['Ped_Inf_R','Ped_Inf_L','Ped_Sup_R','Ped_Sup_L','Plat_Inf_Cent','Plat_Sup_Cent']
    liste_coordonnees = []
    compteur = 0;
    with open (chemin, newline ='') as fichier:
        csvreader = csv.reader(fichier,delimiter=';')
        for row in csvreader: 
            nom = row[0].removeprefix("'").removesuffix("'")
            if nom!= '':
                liste_vertebres.append(row[0].removeprefix("'").removesuffix("'"))    
            valeur_int = row[2].removeprefix("[").removesuffix("]")
            valeur_liste = valeur_int.split(',')
            valeur = [float(i) for i in valeur_liste]
            liste_coordonnees.append(valeur)
    for i in range(len(liste_vertebres)):
        dictionnaire[liste_vertebres[i]] = [liste_coordonnees[i] for i in range(compteur,compteur+6)]
        compteur+=6
    return dictionnaire
#%% Collection des données
Calib_Beads_2Dmat_PA0 = lecture_fichier_2D(r'C:\Users\cleme\Downloads\BDE_Phelma_LOOP_Thème_OLINP_Phelma\Calib_Beads_2Dmat_PA0_csv.csv')
Calib_Beads_2Dmat_PA20 = lecture_fichier_2D(r'C:\Users\cleme\Downloads\BDE_Phelma_LOOP_Thème_OLINP_Phelma\Calib_Beads_2Dmat_PA20_csv.csv')
Calib_Beads_3Dmat = lecture_fichier_2D(r'C:\Users\cleme\Downloads\BDE_Phelma_LOOP_Thème_OLINP_Phelma\Calib_Beads_3Dmat_csv.csv')
Vertebra_2D_PA0 = lecture_fichier_3D(r'C:\Users\cleme\Downloads\BDE_Phelma_LOOP_Thème_OLINP_Phelma\Vertebra_2D_PA0_csv.csv')
Vertebra_2D_PA20 = lecture_fichier_3D(r'C:\Users\cleme\Downloads\BDE_Phelma_LOOP_Thème_OLINP_Phelma\Vertebra_2D_PA20_csv.csv')
#%% liste des 50  billes en communs sur les deux radios 
with open(r'C:\Users\cleme\Downloads\BDE_Phelma_LOOP_Thème_OLINP_Phelma\Name_copy_csv.csv') as fichier:
    csvreader = csv.reader(fichier,delimiter=';')
    for row in csvreader:
        billes_communes = row
#%% Question 1.1
def get_calibration_vector(data_set, billes_3D, billes_a_considerer):
    Mat=[]
    vecteur_u_v = []
    for i in billes_a_considerer:
        u,v=data_set[i][0],data_set[i][1]
        X,Y,Z = billes_3D[i][0],billes_3D[i][1],billes_3D[i][2]
        Mat.append([X,Y,Z,1,0,0,0,0,-u*X,-u*Y,-u*Z])
        Mat.append([0,0,0,0,X,Y,Z,1,-v*X,-v*Y,-v*Z])
        vecteur_u_v.append(u)
        vecteur_u_v.append(v)
    Mat = np.array(Mat)
    Vecteur_uv = np.array(vecteur_u_v)
    ### Calcul du vecteur de calibration
    Vecteur_calibration = np.linalg.inv(np.matrix.transpose(Mat)@Mat)@np.matrix.transpose(Mat)@Vecteur_uv
    # Vecteur_calibration = np.linalg.lstsq(Mat,Vecteur_uv)
    return Vecteur_calibration

Vecteur_calibration_0 = get_calibration_vector(Calib_Beads_2Dmat_PA0,Calib_Beads_3Dmat,billes_communes)
Vecteur_calibration_20 = get_calibration_vector(Calib_Beads_2Dmat_PA20,Calib_Beads_3Dmat,billes_communes)
#%%
#Conversion d'arrays en listes
def array_to_list(un_array):
    liste = []
    for i in range(len(un_array)):
        liste.append(float(un_array[i]))
    return liste
#%%
V_L_PA0 = array_to_list(Vecteur_calibration_0)
V_L_PA20 = array_to_list(Vecteur_calibration_20)
#%% Question 1.2
def get_parameters(Vecteur_calibration):
    dictionnaire ={}
    vec_XYZ = np.array([[Vecteur_calibration[0],Vecteur_calibration[1],Vecteur_calibration[2]],
                        [Vecteur_calibration[4],Vecteur_calibration[5],Vecteur_calibration[6]],
                        [Vecteur_calibration[8],Vecteur_calibration[9],Vecteur_calibration[10]]])@np.array([-Vecteur_calibration[3],Vecteur_calibration[7],-1])
    X0,Y0,Z0= vec_XYZ[0],vec_XYZ[1],vec_XYZ[2]
    d = -1/(np.sqrt((Vecteur_calibration[8]**2)+(Vecteur_calibration[9]**2)+(Vecteur_calibration[10]**2)))
    u0 = (Vecteur_calibration[0]*Vecteur_calibration[8]+Vecteur_calibration[1]*Vecteur_calibration[9]+Vecteur_calibration[2]*Vecteur_calibration[10])*d*d
    v0 = (Vecteur_calibration[4]*Vecteur_calibration[8]+Vecteur_calibration[5]*Vecteur_calibration[9]+Vecteur_calibration[6]*Vecteur_calibration[10])*d*d
    Cu = np.sqrt(d*d*(((u0*Vecteur_calibration[8]-Vecteur_calibration[0])**2)+((u0*Vecteur_calibration[9]-Vecteur_calibration[1])**2)+((u0*Vecteur_calibration[10]-Vecteur_calibration[2])**2)))
    Cv = np.sqrt(d*d*(((v0*Vecteur_calibration[8]-Vecteur_calibration[4])**2)+((v0*Vecteur_calibration[9]-Vecteur_calibration[5])**2)+((v0*Vecteur_calibration[10]-Vecteur_calibration[6])**2)))
    R1_1 = (d/Cu)*(u0*Vecteur_calibration[8]-Vecteur_calibration[0])
    R1_2 = (d/Cu)*(u0*Vecteur_calibration[9]-Vecteur_calibration[1])
    R1_3 = (d/Cu)*(u0*Vecteur_calibration[10]-Vecteur_calibration[2])
    R2_1 = (d/Cv)*(v0*Vecteur_calibration[8]-Vecteur_calibration[4])
    R2_2 = (d/Cv)*(v0*Vecteur_calibration[9]-Vecteur_calibration[5])
    R2_3 = (d/Cv)*(v0*Vecteur_calibration[10]-Vecteur_calibration[6])
    R3_1 = Vecteur_calibration[8]*d
    R3_2 = Vecteur_calibration[9]*d
    R3_3 = Vecteur_calibration[10]*d
    M = [[R1_1,R1_2,R1_3],[R2_1,R2_2,R2_3],[R3_1,R3_2,R3_3]]
    liste_noms = ['X0','Y0','Z0','d','u0','v0','Cu','Cv','Matrice de rotation']
    liste_varaibles = [X0,Y0,Z0,d,u0,v0,Cu,Cv,M]
    for i in range(len(liste_noms)):
        dictionnaire[liste_noms[i]] = liste_varaibles[i]
    return dictionnaire
param_PA0 = get_parameters(V_L_PA0)
parame_PA20 = get_parameters(V_L_PA20)
#%% Question 1.3 reconstruct
def reconstruct(L_PA0,L_PA20,vertebres_PA0,vertebres_PA20):
    liste_X = []
    liste_Y = []
    liste_Z = []
    for V in vertebres_PA0:
        for i in range(6):
            u0,v0 = vertebres_PA0[V][i][0],vertebres_PA0[V][i][1]
            u20,v20 = vertebres_PA20[V][i][0],vertebres_PA20[V][i][1]
            Matrice = np.array([[L_PA0[0]-(L_PA0[8]*u0),L_PA0[1]-(L_PA0[9]*u0),L_PA0[2]-(L_PA0[10]*u0)],
                                [L_PA0[4]-(L_PA0[8]*v0),L_PA0[5]-(L_PA0[9]*v0),L_PA0[6]-(L_PA0[10]*v0)],
                                [L_PA20[0]-(L_PA20[8]*u20),L_PA20[1]-(L_PA20[9]*u20),L_PA20[2]-(L_PA20[10]*u20)],
                                [L_PA20[4]-(L_PA20[8]*v20),L_PA20[5]-(L_PA20[9]*v20),L_PA20[6]-(L_PA20[10]*v20)]])
            Vecteur_c=np.array([u0-L_PA0[3],v0-L_PA0[7],u20-L_PA20[3],v20-L_PA20[7]])
            V_XYZ = np.linalg.inv(np.matrix.transpose(Matrice)@Matrice)@np.matrix.transpose(Matrice)@Vecteur_c
            liste_X.append(V_XYZ[0])
            liste_Y.append(V_XYZ[1])
            liste_Z.append(V_XYZ[2])
    return(liste_X,liste_Y,liste_Z)

X,Y,Z = reconstruct(V_L_PA0,V_L_PA20,Vertebra_2D_PA0,Vertebra_2D_PA20)
ax = plt.figure().add_subplot(projection='3d')
ax.plot([Calib_Beads_3Dmat[i][0] for i in Calib_Beads_3Dmat],[Calib_Beads_3Dmat[i][1] for i in Calib_Beads_3Dmat],[Calib_Beads_3Dmat[i][2] for i in Calib_Beads_3Dmat],'or',)
ax.plot(X,Y,Z,'o')
ax.set_xlabel('X',fontsize=16)
ax.set_ylabel('Y',fontsize=16)
ax.set_zlabel('Z',fontsize=16)    
ax.set_xlim([-600, 100])
ax.set_ylim([-300, 200])
ax.set_zlim([-300, 300])

#%% Question 2.1
def compute_RMS(data_ref,data_test):
    RMS = 0
    taille = len(data_ref)
    for i in range(taille):
        RMS += (data_test[i]-data_ref[i])**2/taille
    RMS_f = np.sqrt(RMS)
    return RMS_f
#print(compute_RMS([1,2,3,4,5,6],[1,2.1,3.1,4.1,5.1,6.1]))

nb_points = []
moyenne_RMS_X = []
moyenne_RMS_Y = []
moyenne_RMS_Z = []
for i in tqdm(range(6,51)):
    nb_point_moy = 200
    RMS_X_int = 0
    RMS_Y_int = 0
    RMS_Z_int = 0
    for j in range(nb_point_moy):
        # liste_points = [random.randint(1,49) for _ in range(i)]
        liste_points = random.sample(list(range(0,50)),i)
        billes_utiles = [billes_communes[k] for k in liste_points]
        L_PA0 = get_calibration_vector(Calib_Beads_2Dmat_PA0,Calib_Beads_3Dmat,billes_utiles)
        L_PA20 = get_calibration_vector(Calib_Beads_2Dmat_PA20,Calib_Beads_3Dmat,billes_utiles)
        l_X,l_Y,l_Z = reconstruct(L_PA0,L_PA20,Vertebra_2D_PA0,Vertebra_2D_PA20)
        l_X_R,l_Y_R,l_Z_R = compute_RMS(X,l_X),compute_RMS(Y,l_Y),compute_RMS(Z,l_Z)
        RMS_X_int += l_X_R
        RMS_Y_int += l_Y_R
        RMS_Z_int += l_Z_R
    nb_points.append(i)
    moyenne_RMS_X.append(RMS_X_int/nb_point_moy)
    moyenne_RMS_Y.append(RMS_Y_int/nb_point_moy)
    moyenne_RMS_Z.append(RMS_Z_int/nb_point_moy)
#%% Tracé des erreurs
plt.plot(nb_points,moyenne_RMS_X, label = 'RMS X')
plt.plot(nb_points,moyenne_RMS_Y, label = 'RMS Y')
plt.plot(nb_points,moyenne_RMS_Z, label = 'RMS Z')
plt.legend(fontsize = 18)
plt.style.use('seaborn')
plt.xlabel("Nombre de billes utilisées pour la reconstruction",fontsize=18)
plt.ylabel("Erreur RMS",fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title("Evolution de l'erreur RMS en fonction du nombre de points de calibration")
plt.grid()
plt.show()

#%% Question 2.2
def get_volume(indices_A,indices_B):
    coord_A = [billes_communes[i] for i in indices_A]
    PlanA = [Calib_Beads_3Dmat[i] for i in coord_A]
    coord_B = [billes_communes[i] for i in indices_B]
    PlanB = [Calib_Beads_3Dmat[i] for i in coord_B]
    surface_1 = 0.5*np.abs((PlanA[0][1]*(PlanA[1][2]-PlanA[2][2])+PlanA[1][1]*(PlanA[2][2]-PlanA[0][2])+PlanA[2][1]*(PlanA[0][2]-PlanA[1][2])+PlanA[3][1]*(PlanA[1][2]-PlanA[0][2])))
    surface_2 = 0.5*np.abs((PlanB[0][1]*(PlanB[1][2]-PlanB[2][2])+PlanB[1][1]*(PlanB[2][2]-PlanB[0][2])+PlanB[2][1]*(PlanB[0][2]-PlanB[1][2])+PlanB[3][1]*(PlanB[1][2]-PlanB[0][2])))
    return((surface_1 + surface_2)*300)
#%%
# Dans la liste "billes communes", les billes A correspondent occupent les indices 0 à 25
Volume = []
RMS_X_vol = []
RMS_Y_vol = []
RMS_Z_vol = []
for j in tqdm(range(20000)):
    liste_pointsA = random.sample(list(range(0,26)),4)
    liste_pointsB = random.sample(list(range(0,50)),4)
    billes_utiles_vol = []
    billes_a = []
    billes_b = []
    for i in range(4):
        billes_a.append(liste_pointsA[i])
        billes_b.append(liste_pointsB[i])
    i_volume = get_volume(billes_a,billes_b)
    billes_utiles_vol = billes_a+billes_b
    billes_utiles_vf = [billes_communes[k] for k in billes_utiles_vol]
    L_PA0 = get_calibration_vector(Calib_Beads_2Dmat_PA0,Calib_Beads_3Dmat,billes_utiles_vf)
    L_PA20 = get_calibration_vector(Calib_Beads_2Dmat_PA20,Calib_Beads_3Dmat,billes_utiles_vf)
    l_X,l_Y,l_Z = reconstruct(L_PA0,L_PA20,Vertebra_2D_PA0,Vertebra_2D_PA20)
    l_X_R,l_Y_R,l_Z_R = compute_RMS(X,l_X),compute_RMS(Y,l_Y),compute_RMS(Z,l_Z)
    Volume.append(i_volume)
    RMS_X_vol.append(l_X_R)
    RMS_Y_vol.append(l_Y_R)
    RMS_Z_vol.append(l_Z_R)
#%%
# Maintenant on trie les listes selon les volumes décroissants
zipped_lists = list(zip(Volume,RMS_X_vol,RMS_Y_vol,RMS_Z_vol))
zipped_lists.sort(key=lambda x: x[0])
liste_vol_triee, liste_RMSX_triee, liste_RMSY_triee, liste_RMSZ_triee = zip(*zipped_lists)
plt.plot(liste_vol_triee,liste_RMSX_triee,label='RMS X')
plt.plot(liste_vol_triee,liste_RMSY_triee,label='RMS Y')
plt.plot(liste_vol_triee,liste_RMSZ_triee,label='RMS Z')
plt.legend(fontsize = 18)
plt.style.use('seaborn')
plt.xlabel("Volume",fontsize=18)
plt.ylabel("Erreur RMS",fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title("Evolution de l'erreur RMS en fontion du volume",fontsize=18)
plt.grid()
plt.show()
#%% Question 2.3
def compute_euclidien_distance(pointA,pointB):
    d = np.sqrt(((pointA[0]-pointB[0])**2)+((pointA[1]-pointB[1])**2)+((pointA[2]-pointB[2])**2))
    return d 
# print(compute_euclidien_distance([3,4,5],[6,5,8]))
#%%
billes_2_3 = ['A_2_3', 'A_2_5', 'A_3_3', 'A_3_5', 'B_2_2', 'B_2_3', 'B_3_2', 'B_3_3']
calibration_volume = []
for i in billes_2_3:
    calibration_volume.append(Calib_Beads_3Dmat[i])
barycentre_volume = []
for i in range(3):
    barycentre_volume.append(np.mean([calibration_volume[k][i] for k in range(len(calibration_volume))]))
Vecteur_calibration_0_2_3 = get_calibration_vector(Calib_Beads_2Dmat_PA0,Calib_Beads_3Dmat,billes_2_3)
Vecteur_calibration_20_2_3 = get_calibration_vector(Calib_Beads_2Dmat_PA20,Calib_Beads_3Dmat,billes_2_3)
X_2_3,Y_2_3,Z_2_3 = reconstruct(Vecteur_calibration_0_2_3,Vecteur_calibration_20_2_3,Vertebra_2D_PA0,Vertebra_2D_PA20)
distance_vertebres = [compute_euclidien_distance([X[i],Y[i],Z[i]],[X_2_3[i],Y_2_3[i],Z_2_3[i]]) for i in range(len(X))]
distance_barycentre = [compute_euclidien_distance([barycentre_volume[0],barycentre_volume[1],barycentre_volume[2]],[X_2_3[i],Y_2_3[i],Z_2_3[i]]) for i in range(len(X))]
#%%
plt.plot(distance_barycentre,distance_vertebres,'o')
plt.legend(fontsize = 18)
plt.style.use('seaborn')
plt.xlabel("Distance par rapport au centre de gravité",fontsize=18)
plt.ylabel("Distance entre deux points corresndants",fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title("Evolution de l'erreur de reconstruction en fonction de la distance par rapport au barycentre",fontsize=18)
plt.grid()
plt.show()
#%% ajout du bruit Gaussien

def noise(sigma, dict_a_bruiter,dim):

    # On initialise des listes contenant les coordonnées u et v de chaque bead
    dict_bruite=copy.deepcopy(dict_a_bruiter)
    
    if dim==2:
        u_beads=[]
        v_beads=[]    
        # On acquiert ces coordonnées
        for k in row:
            u_beads.append(dict_a_bruiter[k][0])
            v_beads.append(dict_a_bruiter[k][1])
            
        for k in row:
            dict_bruite[k][0] = random.gauss(dict_bruite[k][0],sigma)
            dict_bruite[k][1] = random.gauss(dict_bruite[k][1],sigma)
            
        
    
    elif dim==3:
        x_beads=[]
        y_beads=[] 
        z_beads=[]
        # On acquiert ces coordonnées
        for k in row:
            x_beads.append(dict_a_bruiter[k][0])
            y_beads.append(dict_a_bruiter[k][1])
            z_beads.append(dict_a_bruiter[k][2])

        for k in row:
            dict_bruite[k][0] = random.gauss(dict_bruite[k][0],sigma)
            dict_bruite[k][1] = random.gauss(dict_bruite[k][1],sigma)
            dict_bruite[k][2] = random.gauss(dict_bruite[k][2],sigma)


    
    return dict_bruite

#%% Question 3.1
sigma=[0.05,0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15]
longueur_sig=len(sigma)
nb_tot_rep=200
l_x_r_noisy_arr= np.zeros((longueur_sig,nb_tot_rep), dtype = float)
l_y_r_noisy_arr= np.zeros((longueur_sig,nb_tot_rep), dtype = float)
l_z_r_noisy_arr= np.zeros((longueur_sig,nb_tot_rep), dtype = float)

c=0;
for rep in range(nb_tot_rep):
    for j in sigma:
        Calib_Beads_2Dmat_PA0_noisy=noise(j,Calib_Beads_2Dmat_PA0,2)
        Calib_Beads_2Dmat_PA20_noisy=noise(j,Calib_Beads_2Dmat_PA20,2)
     #   Calib_Beads_3Dmat_noisy=noise(j,Calib_Beads_3Dmat,3)
        
        
        L_PA0_noisy = get_calibration_vector(Calib_Beads_2Dmat_PA0_noisy,Calib_Beads_3Dmat,billes_utiles)
        L_PA20_noisy = get_calibration_vector(Calib_Beads_2Dmat_PA20_noisy,Calib_Beads_3Dmat,billes_utiles)
        l_X_noisy,l_Y_noisy,l_Z_noisy = reconstruct(L_PA0_noisy,L_PA20_noisy,Vertebra_2D_PA0,Vertebra_2D_PA20)
        # ax = plt.figure().add_subplot(projection='3d')
        # ax.plot(X,Y,Z,'ob')
        # ax.plot(l_X_noisy,l_Y_noisy,l_Z_noisy, 'or')
        # ax.set_xlabel('X',fontsize=16)
        # ax.set_ylabel('Y',fontsize=16)
        # ax.set_zlabel('Z',fontsize=16)    
        # ax.set_xlim([-600, 100])
        # ax.set_ylim([-300, 200])
        # ax.set_zlim([-300, 300])
        
        l_x_r_noisy_arr[sigma.index(j),rep],l_y_r_noisy_arr[sigma.index(j),rep],l_z_r_noisy_arr[sigma.index(j),rep] = compute_RMS(X,l_X_noisy),compute_RMS(Y,l_Y_noisy),compute_RMS(Z,l_Z_noisy)
    

list_err_x=[]
list_err_y=[]
list_err_z=[]
for k in range(len(sigma)):
    list_err_x.append(np.mean(l_x_r_noisy_arr[k,:]))
    list_err_y.append(np.mean(l_y_r_noisy_arr[k,:]))
    list_err_z.append(np.mean(l_z_r_noisy_arr[k,:]))
    

plt.plot(sigma,list_err_x, label='RMS X')
plt.plot(sigma,list_err_y, label='RMS Y')
plt.plot(sigma,list_err_z, label='RMS Z')
plt.legend(fontsize = 18)
plt.xlabel("Valeur de sigma",fontsize=18)
plt.ylabel("RMSE",fontsize=18)
plt.title("Evolution de l'erreur de reconstruction en fonction de la variance du bruit ajouté",fontsize=18)

#%% Question 3.2

sigma=[0.05,0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15]
longueur_sig=len(sigma)
nb_tot_rep=200
l_x_r_noisy_arr= np.zeros((longueur_sig,nb_tot_rep), dtype = float)
l_y_r_noisy_arr= np.zeros((longueur_sig,nb_tot_rep), dtype = float)
l_z_r_noisy_arr= np.zeros((longueur_sig,nb_tot_rep), dtype = float)

c=0;
for rep in range(nb_tot_rep):
    for j in sigma:
        Calib_Beads_2Dmat_PA0_noisy=noise(j,Calib_Beads_2Dmat_PA0,2)
        Calib_Beads_2Dmat_PA20_noisy=noise(j,Calib_Beads_2Dmat_PA20,2)
     #   Calib_Beads_3Dmat_noisy=noise(j,Calib_Beads_3Dmat,3)
        
     
        liste_pointsA = random.sample(list(range(0,26)),4)
        liste_pointsB = random.sample(list(range(0,50)),4)
        billes_utiles_vol = []
        billes_a = []
        billes_b = []
        for i in range(4):
            billes_a.append(liste_pointsA[i])
            billes_b.append(liste_pointsB[i])
        i_volume = get_volume(billes_a,billes_b)
        billes_utiles_vol = billes_a+billes_b
        billes_utiles_vf = [billes_communes[k] for k in billes_utiles_vol]
        
        L_PA0_noisy = get_calibration_vector(Calib_Beads_2Dmat_PA0_noisy,Calib_Beads_3Dmat,billes_utiles_vf)
        L_PA20_noisy = get_calibration_vector(Calib_Beads_2Dmat_PA20_noisy,Calib_Beads_3Dmat,billes_utiles_vf)
        l_X_noisy,l_Y_noisy,l_Z_noisy = reconstruct(L_PA0_noisy,L_PA20_noisy,Vertebra_2D_PA0,Vertebra_2D_PA20)
        # ax = plt.figure().add_subplot(projection='3d')
        # ax.plot(X,Y,Z,'ob')
        # ax.plot(l_X_noisy,l_Y_noisy,l_Z_noisy, 'or')
        # ax.set_xlabel('X',fontsize=16)
        # ax.set_ylabel('Y',fontsize=16)
        # ax.set_zlabel('Z',fontsize=16)    
        # ax.set_xlim([-600, 100])
        # ax.set_ylim([-300, 200])
        # ax.set_zlim([-300, 300])
        
        l_x_r_noisy_arr[sigma.index(j),rep],l_y_r_noisy_arr[sigma.index(j),rep],l_z_r_noisy_arr[sigma.index(j),rep] = compute_RMS(X,l_X_noisy),compute_RMS(Y,l_Y_noisy),compute_RMS(Z,l_Z_noisy)
    

list_err_x=[]
list_err_y=[]
list_err_z=[]
for k in range(len(sigma)):
    list_err_x.append(np.mean(l_x_r_noisy_arr[k,:]))
    list_err_y.append(np.mean(l_y_r_noisy_arr[k,:]))
    list_err_z.append(np.mean(l_z_r_noisy_arr[k,:]))
    

plt.plot(sigma,list_err_x, label='RMS X')
plt.plot(sigma,list_err_y, label='RMS Y')
plt.plot(sigma,list_err_z, label='RMS Z')
plt.legend(fontsize = 18)
plt.xlabel("Valeur de sigma",fontsize=18)
plt.ylabel("RMSE",fontsize=18)
plt.title("Evolution de l'erreur de reconstruction en fonction de la variance du bruit ajouté",fontsize=18)
