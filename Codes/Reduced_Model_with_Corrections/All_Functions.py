# -*- coding: utf-8 -*-
"""
FUNCTIONS
"""

from numpy import *
#import numpy as np
from decimal import Decimal
import numpy as np











#################### GRID_COORDINATES #########################################
"""
INPUTS:
- nastran_file: a bdf file containing grid points (card GRID)
- list_nodes: a list of nodes corresponding to some of those grid points

OUTPUTS:
- coord_grid: a dictionary with each node in the list and its coordinates:
{node : [x_coord,y_coord,z_coord]}
"""
def grid_coordinates(nastran_file, list_nodes): 
    # Nastran file reading
    datfile = open(nastran_file,'r')
    text    = datfile.read()
    lines   = text.splitlines()
    datfile.close()
    # Creation of the dictionary with nodes and coordinates 
    coord_grid = {}
    for line in lines:
        if 'GRID' in line:
            if ',' in line:
                line_splitted = line.split(',')
            else:
                line_splitted = line.split()
            grid_node = int(line_splitted[1])
            if grid_node in list_nodes:
                grid_coordinates      = [line_splitted[3],line_splitted[4],line_splitted[5]]
                coord_grid[grid_node] = map(float, grid_coordinates)
            
    return coord_grid





######################## FREE #################################################
# function to format card entries to fit in the 8 digits of free field nastran format
# output is a string with a number occupying 8 characters
def free(number):
    # Number of characters of the number
    n_characters = len(str(number))
    # The formating is applied only if the number takes more than 8 characters
    if n_characters > 8:    
        # Exponent of the number (e.g. 1325 -> exponent = 3)    
        exponent = np.floor(np.log10(np.abs(number)))
        # If the number is negative, one character will be the minus sign
        # 8 characters:  exp<10 -> -2.456+8  , exp>10 ->  -2.45+78
        if number < 0:
            if abs(exponent) < 10:
                number_free_format = "{:.3E}".format(Decimal(str(number)))
            else:
                number_free_format = "{:.2E}".format(Decimal(str(number)))
        # If the number is positive, there s no minus sign
        # 8 characters:  exp<10 -> 1.3456+8  , exp>10 ->  1.345+78
        else:
            if abs(exponent) < 10:
                number_free_format = "{:.4E}".format(Decimal(str(number)))
            else:
                number_free_format = "{:.3E}".format(Decimal(str(number)))        

    else:
        number_free_format = str(number)
        
    # Finally character E is deleted as it is not needed in nastran 
    number_free_format = number_free_format.replace("E", "")
    number_free_format = number_free_format.replace("e", "") # NEW
        
    return number_free_format










######################## WriteForceMoment #####################################
# (from Osycaf project)

##########################
# \brief Écriture de cartes FORCE et MOMENT dans un .dat Nastran.
# \details Écrit les valeurs des cartes FORCE et MOMENT dans un fichier Nastran à partir d'un tableau de données (array).
# \param F un \e array contenant les valeurs des forces et moments ainsi que les numéros des noeuds d'application (voir note pour le format).
# \param [coord_system] le numéro du repère d'écriture des forces. Par défaut 0.
# \param [ddl] la liste du choix des ddl ('TX','TY','TZ','RX','RY','RZ') à écrire. La valeur est nulle pour les ddl non inclus dans la liste.
# \param [common_load] un booléen pour savoir si toutes les cartes doivent être regroupées dans une seule carte LOAD. Par défaut True. Dnas le cas False, une carte LOAD est créée pour chaque carte FORCE et MOMENT et la liste des identifiants des cartes LOAD est retournée en sortie.
# \return La chaine de caractère contenant toutes les cartes, à inclure directement dans le fichier .dat Nastran.
# \return L'identifiant de la carte LOAD (ou la liste des identifiants si common_load=False).
# \note \li La carte LOAD (combinaison linéaire des efforts et moments) est automatiquement créée avec des coefficients unitaires. Les identifiants des cartes FORCE et MOMENT sont numérotés à partir de 1. L'identifiant de la carte LOAD vaut l'identifiant max +1. Il faut penser à ajuster l'identifiant de la carte LOAD dans le CASE CONTROL. \n
# \li Le tableau \e F doit posséder le data type suivant : \code
# dtype([('id',int,1),('value',float64,(6,))]) \endcode pour les 6 ddl dans l'ordre.

def WriteForceMoment(F, coord_system=0, ddl=['TX','TY','TZ','RX','RY','RZ'], common_load=True):

	forces_card = ''	# chaine de caractère vide pour les cartes FORCE
	compt_id = 1	# initialisation du compteur des identifiants 	

	# liste binaire pour le choix des ddl
	f=[]
	m=[]
	ddlf=['TX','TY','TZ']
	ddlm=['RX','RY','RZ']
	for i in range(3):
		if ddlf[i] in ddl:
			f.append(1.)
		else:
			f.append(0.)
		if ddlm[i] in ddl:
			m.append(1.)
		else:
			m.append(0.)
	f=array(f)
	m=array(m)

	for line in F:
		# FORCE
		# dans la cas où les valeurs X, Y et Z sont nulles (pb nastran)
		array_value=line['value'][0:3]*f
		if (array_value == zeros(3)).all():
			F_value = 0.
			array_value = ones(3)
		else:
			F_value=1.
		forces_card += 'FORCE   {0:8}{1:8}{2:8}{3:8}{4:8}{5:8}{6:8}\n'.format(*[int(compt_id), int(line['id']), \
										int(coord_system), F_value, array_value[0], \
										array_value[1], array_value[2]])

		# MOMENT
		# dans la cas où les valeurs X, Y et Z sont nulles (pb nastran)
		array_value=line['value'][3:6]*m
		if (array_value == zeros(3)).all():
			M_value = 0.
			array_value = ones(3)
		else:
			M_value=1.

		forces_card += 'MOMENT  {0:8}{1:8}{2:8}{3:8}{4:8}{5:8}{6:8}\n'.format(*[int(compt_id), int(line['id']), \
										int(coord_system), M_value, array_value[0], \
										array_value[1], array_value[2]])
		compt_id += 1

	forces_card += '\n\n'
	
	# ecriture carte LOAD
	if common_load==True:
		total_id = compt_id
		compt_id = 1
		forces_card += 'LOAD    {0:8}{1:8}{2:8}{3:8}{4:8}{5:8}{6:8}{7:8}\n'.format(*[total_id+1, 1., 1., \
												compt_id, 1., compt_id+1, 1., \
												compt_id+2])
		compt_id +=3
		for i in range(compt_id,total_id,4):
			forces_card += '        {0:8}{1:8}{2:8}{3:8}{4:8}{5:8}{6:8}{7:8}\n'.format(*[1., i, 1., i+1, 1., i+2, 1., i+3])
		id_load = total_id+1
	else:
		nb_id=compt_id
		id_load=[]
		for i in range(1,nb_id):
			compt_id+=1
			forces_card += 'LOAD    {0:8}{1:8}{2:8}{3:8}\n'.format(*[compt_id, 1., 1., i])
			id_load.append(compt_id)

	return (forces_card, id_load)




