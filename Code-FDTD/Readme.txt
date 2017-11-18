1/ Avant le calcul, modifier les chemins de l'ordinateur dans les programmes 
'FDTD_3D_ripples_LD_main' et "load_E_last_half_an_optical_cycle_everytimestep" et "load_data_z"
Le chemin 'exci' servira de stockage des données calculées

2/ Lancer d'abord  "FDTD_3D_ripples_LD_main"

3/ Lancer ensuite "load_E_last_half_an_optical_cycle_everytimestep" et cliquer pour faire défiler le champ

4/ Enfin "load_data_z" permet l'affichage de l'énergie absorbée (entre autres) 

5/ Attention, dans le calcul le champ est polarisé suivant "z" et la propagation se fait selon "x". Dans l'affichage l'axe des absisses correspond à z mais est appelé x (les deux axes sont inversés uniquement pour l'affichage) 

6/ pulse_number=1; Tester 1 ou 2 pulses

7/ Augmenter la dimension de la boîte de calcul en fonction de la RAM de l'ordinateur
Y_dimension_main=12e-6; 
Z_dimension_main=12e-6;

8/ roughness_type=1; Créer un type 3 pour une rugosité bien définie, calculée dans un programme externe "Roughness_defined"

9/ La rugosité est actuellement définie dans la section "roughness_function" du programme "FDTD_3D_ripples_LD_loop"

10/ Faire une étude comparative de l'effet de la rugosité sur la FFT de l'énergie absorbée (via "load_data_z") et étudier l'évolution du champ électromagnétique associé (via "load_E_last_half_an_optical_cycle_everytimestep") .
