1/ Avant le calcul, modifier les chemins de l'ordinateur dans les programmes 
'FDTD_3D_ripples_LD_main' et "load_E_last_half_an_optical_cycle_everytimestep" et "load_data_z"
Le chemin 'exci' servira de stockage des donn�es calcul�es

2/ Lancer d'abord  "FDTD_3D_ripples_LD_main"

3/ Lancer ensuite "load_E_last_half_an_optical_cycle_everytimestep" et cliquer pour faire d�filer le champ

4/ Enfin "load_data_z" permet l'affichage de l'�nergie absorb�e (entre autres) 

5/ Attention, dans le calcul le champ est polaris� suivant "z" et la propagation se fait selon "x". Dans l'affichage l'axe des absisses correspond � z mais est appel� x (les deux axes sont invers�s uniquement pour l'affichage) 

6/ pulse_number=1; Tester 1 ou 2 pulses

7/ Augmenter la dimension de la bo�te de calcul en fonction de la RAM de l'ordinateur
Y_dimension_main=12e-6; 
Z_dimension_main=12e-6;

8/ roughness_type=1; Cr�er un type 3 pour une rugosit� bien d�finie, calcul�e dans un programme externe "Roughness_defined"

9/ La rugosit� est actuellement d�finie dans la section "roughness_function" du programme "FDTD_3D_ripples_LD_loop"

10/ Faire une �tude comparative de l'effet de la rugosit� sur la FFT de l'�nergie absorb�e (via "load_data_z") et �tudier l'�volution du champ �lectromagn�tique associ� (via "load_E_last_half_an_optical_cycle_everytimestep") .
