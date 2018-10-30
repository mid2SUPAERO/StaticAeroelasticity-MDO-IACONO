# -*- coding: utf-8 -*-



#from numpy import *
import All_Functions as function
#from numpy import *
import numpy as np














#################### GENERATE_LOADS ###########################################
"""
DESCRIPTION =================================================================

INPUTS
- A list of nodes, from clamped end to tip

- A nastran input file (.bdf) which should contain
Two pairs of markers $${ $$} to specify where to write things on the bdf
In the 1st pair it will write SET and SUBCASE cards
In the 2nd pair it will write FORCE,MOMENT,LOAD cards

- At the 1st node (which is clamped) loads are not applied
but displacements are measured


OUTPUTS
- It will write the cards for 6 unit loads at each node in the list 
(except at the 1st one)
- It will write as output request displacements at list of nodes
"""

# What this program will do is to generate different load cases (subcases):
# for each 'central', 6 load cases will be generated,
# (a unit effort in each component Fx,Fy,Fz,Mx,My,Mz)
# To generate a load case, the following cards are written on the bdf:
# SUBCASE,FORCE,MOMENT,LOAD
# When running 'nastran_file_mod' a punch file will be created
# It will contain the displacement at the 'centrals' for each subcase

# Note:
# A single bdf contains all the analyses needed
# A single punch file contains all results


## List of 'centrals' 
#list_nodes = [7000,7001,7002,7003,7004,7005,7006,7007,7008,7009,7010]
#
## bdf with the characteristics described above
#nastran_file = "goland_wing_extraction.bdf"
#
## bdf where the modified bdf will be written
#nastran_file_mod = "goland_wing_extraction_mod.bdf"





def generate_loads(list_nodes,nastran_file,S_Ref):

    nastran_file_mod = nastran_file[:-4]+'M.inp'
    
    # Careful! in python if you do a=b
    # modifications in b will affect a (and viceversa)
    # to do a copy you must do a=b[:]
    
    # 1st node is conserved for the set where displacements are measured
    nodes_set = list_nodes[:]
    
    # And it s deleted for node where loads are applied
    list_nodes.pop(0)
    
    #ddl=['TX','TY','TZ','RX','RY','RZ']
    
    
    # Description of the structure of Nastran cards
    # FORCE  - Load ID - Node of application - Coord. System - Scale factor - Fx - Fy - Fz   
    # MOMENT - Load ID - Node of application - Coord. System - Scale factor - Mx - My - Mz 
    # LOAD - Load Set ID - Overall scale factor - Pairs of (Load scale factor, Load ID)
    # Note: FORCE and MOMENT ID can be the same
    # Note: LOAD ID must be different from the FORCE and MOMENT IDs
    
    
    # Creation of unitary efforts (forces and moments) table ======================
    # It creates a list whose elements have the following structure:
    # (node,[Fx,Fy,Fz,Mx,My,Mz])
    # This will be used by function.WriteForceMoment function
    
    forces=[] 
    
    for node in list_nodes: 
        for k in range(6): #create six cases for each node (one unitary effort, the rest null)
            f=np.zeros(6)
            f[k]=1.
            forces.append((int(node), f))
    forces=np.array(forces, dtype=np.dtype([('id',int),('value',float,(6,))])) #define the data types (node=integer,effort componenets=float)
    
    
    # Creation of nastran FORCE,MOMENT and LOAD cards =============================
    forces_card, list_id_load=function.WriteForceMoment(forces, coord_system=S_Ref, \
            		ddl=['TX','TY','TZ','RX','RY','RZ'], common_load=False)
    
    
    # Creation of nastran SET card ================================================
    # Creation of the SET containing the 'centrals'
    set_node='SET 1='
    cl=0
    for n in nodes_set:
    	if cl==10:
    		set_node+='\n'
    		cl=0
    	set_node+=str(n)+','
    	cl+=1		
    set_node=set_node[0:-1]+'\n\n'
    
    # Creation of nastran SUBCASE cards ===========================================        
    subcases=''
    subcases+='\nDISP(PUNCH)=1\n'
    for i,sub in enumerate(list_id_load):
    	subcases+='SUBCASE '+str(i+1)+'\n'
    	subcases+='    LOAD='+str(sub)+'\n'
    
    # Writing the created cards on the bdf ========================================
    
    # Nastran file reading
    datfile=open(nastran_file,'r')
    text = datfile.read()
    datfile.close()
    # We save the whole bdf in variable 'text'
    # We edit 'text' and then we write it in the bdf
    
    # Writing SET,SUBCASE cards in CASE CONTROL
    start_first_tag_index = text.find('$${')	# position tag start of CASE CONTROL
    end_first_tag_index   = text.find('$$}')	# position tag  end  of CASE CONTROL
    text = text.replace(text[start_first_tag_index: end_first_tag_index+3], '$${\n'+set_node+subcases+'$$}',1)	
    
    # Writing FORCE,MOMENT,LOAD cards in BULK
    start_second_tag_index = text.rfind('$${')	# position tag start of BULK
    end_second_tag_index   = text.rfind('$$}')	# position tag  end  of BULK
    text = text.replace(text[start_second_tag_index: end_second_tag_index+3], '$${\n'+forces_card+'$$}')
    
    # Writing on the file and closing
    datfile=open(nastran_file_mod,'w')
    datfile.write(text)
    datfile.close()


    return nastran_file_mod


























#################### READ_PUNCH ###############################################
"""
===== DESCRIPTION =====
Code reading deflections (displacements and rotations) of different
subcases in a punch file
It reads them and place them in a dictionary

Here implemented as a function called: deflection_read



===== DETAILED DESCRIPTION =====


===== INPUT =====
name_punch -> a punch file 

The punch file contains the following data structure repeated for each subcase

$TITLE   = MSC.NASTRAN JOB CREATED ON 20-APR-17 AT 15:25:24                  343
$SUBTITLE=                                                                   344
$LABEL   =                                                                   345
$DISPLACEMENTS                                                               346
$REAL OUTPUT                                                                 347
$SUBCASE ID =          20                                                    348
      1021       G      0.000000E+00      0.000000E+00      0.000000E+00     349
-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00     350
      1022       G      1.887874E-03      2.545759E+00     -7.569344E-05     351
-CONT-                 -7.973458E-03      4.664505E-06     -4.299313E-03     352
      1023       G      6.608060E-03      8.972890E+00     -1.212837E-04     353
-CONT-                 -1.280523E-02      7.463462E-06     -1.518801E-02     354
      1024       G     -3.153043E-02      1.588142E+01      2.774670E-02     355
-CONT-                 -1.337688E-02      2.381069E-04     -2.185172E-02     356
      1025       G      1.907242E-02      2.603684E+01     -1.364872E-04     357
-CONT-                 -1.441751E-02      8.316507E-06     -4.403515E-02     358
      1026       G      0.000000E+00      0.000000E+00      0.000000E+00     359
-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00     360




===== OUTPUT =====
A dictionary with the following structure:
note: use this code to access results
#print(DISPLACEMENT)
#print(DISPLACEMENT['SUBCASE1'])
#print(DISPLACEMENT['SUBCASE1'][1025])
#print(DISPLACEMENT['SUBCASE1'][1025]['TX'])

in this example:
'SUBCASE1' -> subcase id
1025       -> node where displacements are measured
'TX'       -> degree of freedom to measure: ['TX','TY','TZ','RX','RY','RZ']

"""


def read_punch(name_punch): 
  
     
    # About text processing in python
    # line=lec_file.readline() => reads the next  line
    # line.split() => it splits the line into words, deleting spaces
    
    # INPUT: punch file to read
    #name_punch = 'check.pch'
    #name_punch = 'beam3d_good.pch'
    #name_punch = 'main_input_mod.pch'
    
    # Obtaining the number of lines in the text (necessary for the for loop below)
    lec_file  = open(name_punch,'r')
    text      = lec_file.read()
    lines     = text.splitlines()
    num_lines = len(lines)
    lec_file.close()
    
    # We open the file again, to read it
    lec_file = open(name_punch,'r')
    
    # line of type 1: contains node number and displacements TX,TY,TZ
    #       1023       G      1.745092E+01     -2.801950E-02      9.427474E-05      11
    # line of type 2: contains  rotations RX,RY,RZ
    # -CONT-                  3.839260E-05      2.489215E-02     -1.996405E-02      12
    
    DISPLACEMENT = {}
    
    for i in range(num_lines):
    
        line = lec_file.readline()
    
        # $SUBCASE indicates the start of the subcase
        # we get the subcase ID
        if '$SUBCASE' in line:
            subcase_id = line.split()[3]
            subcase    =  'SUBCASE'+str(subcase_id)
            DISPLACEMENT[subcase] = {}
            
            # we enter the loop with a line of type 1        
            line=lec_file.readline()
            
            # When the subcase ends, the next line is $TITLE
            # Therefore we keep reading until the word $TITLE appears
            # When we read the last subcase, the following line will be empty
            # instead of containing $TITLE
            while ('$TITLE' not in line) and (line != ''):
    
                # We read the line of type 1
                num_node = int(line.split()[0])
                tx       = float(line.split()[2])
                ty       = float(line.split()[3])
                tz       = float(line.split()[4])
                # We read the line of type 2
                line = lec_file.readline()
                rx   = float(line.split()[1])
                ry   = float(line.split()[2])
                rz   = float(line.split()[3])
                
                # We place the deflections we have read into the dictionary
                DISPLACEMENT[subcase][num_node]       = {}
                DISPLACEMENT[subcase][num_node]["TX"] = tx
                DISPLACEMENT[subcase][num_node]["TY"] = ty
                DISPLACEMENT[subcase][num_node]["TZ"] = tz
                DISPLACEMENT[subcase][num_node]["RX"] = rx
                DISPLACEMENT[subcase][num_node]["RY"] = ry
                DISPLACEMENT[subcase][num_node]["RZ"] = rz
                
                # We move to the next line 
                # (which will be of type 1, contain $TITLE or be empty)
                line = lec_file.readline()
         
    return DISPLACEMENT














#################### BUILD_FLEXIBILITY_MATRICES ###############################
"""
Extraction of section properties of one box

Inputs: 
- punch file 
- Nodes i and j of the box
- Subcases corresponding to 6 units loads applied at node j


Outputs:
- flexibility matrices FF_i, FF_j
Flexibility Matrix = FF
dimension: 6x6
FFi => Relation between force at j and deflections at i
FFj => Relation between force at j and deflections at j
"""

def fleximatrices(name_punch, node_i, node_j, subcases): 
 
    DISPLACEMENT = read_punch(name_punch)
    # Output strucutre: (a dictionary)
    # DISPLACEMENT['SUBCASE1'][1025]['TX']
    
    FF_i = np.zeros((6,6))
    FF_j = np.zeros((6,6))
    dofs = ['TX','TY','TZ','RX','RY','RZ']
   
    for j,subcase in enumerate(subcases):
        for i,dof in enumerate(dofs):
            FF_i[i,j] = DISPLACEMENT[subcase][node_i][dof]
            FF_j[i,j] = DISPLACEMENT[subcase][node_j][dof]
            
    return(FF_i, FF_j)


























#################### EXTRACTION_FUNCTION ######################################



def properties_extraction(D1,D2,l1,l2,E,G):
    
    
    """
    Function extraction
    From a generic box in the full 3D FEM model equivalent beam section properties
    are obtained
    Reference: Cirillo thesis
    
    x as beam axis
    
    Inputs:
    - D1,D2: Flexibility matrices at stations 1 and 2 of the generic box
    These matrices are computed from static analyses in Nastran,
    applying unit forces and moments
    
    - E,G: material properties
    
    - l1,l2: Distances of the box sections from the fixed beam root
    
    
    Outputs:
    - y_NA,z_NA : Axial center coordinates (center of gravity)
    - y_SC,z_SC : Shear center coordinates
    - a         : Section area
    - Iyy,Izz   : Section inertias
    - J         : Section torsional modulus
    - Ky,Kz     : Reduced shear area factor 
    - F         : Flexibility matrix of the box modified (ux,uy only shear deformation)
    - K         : Stiffness matrix of the box modified   (ux,uy only shear deformation)
    """
    
    import numpy as np
    
    
    
    
    # Stiffness matrix Forces - deformations ======================================
    
    # The stiffness matrix should have the following shape
    # |F_x|   |A  0  0  0  A  A|  |epsilon_x|
    # |F_y|   |0  S  S  S  0  0|  |gamma_y  |
    # |F_z| = |0  S  S  S  0  0|· |gamma_z  |
    # |M_x|   |0  S  S  S  0  0|  |kappa_x  |
    # |M_y|   |A  0  0  0  A  A|  |kappa_y  |
    # |M_z|   |A  0  0  0  A  A|  |kappa_z  |
    
    # The zeros in the flexibility matrix should be the same
    # |epsilon_x|   |F  0  0  0  F  F|  |F_x|
    # |gamma_y  |   |0  F  F  F  0  0|  |F_y|
    # |gamma_z  | = |0  F  F  F  0  0|· |F_z|
    # |kappa_x  |   |0  F  F  F  0  0|  |M_x|
    # |kappa_y  |   |F  0  0  0  F  F|  |M_y|
    # |kappa_z  |   |F  0  0  0  F  F|  |M_z|
    
    
    # Initial Flexibiblity Matrix of the box
    F = (D2-D1)/(l2-l1)
    
    # Some elements are forced to be zero
    # (check if they are close to 0 in a real case)
    # Python code    # Mathematical notation
    F[ 0  , 1:4] = 0 #row  1,    , columns 2,3,4
    F[1:4 ,  0 ] = 0 #rows 2,3,4 , column  1
    F[1:4 , 4:6] = 0 #rows 2,3,4 , columns 5,6
    F[4:6 , 1:4] = 0 #rows 5,6   , columns 2,3,4
    
    
    # Initial Stiffness Matrix of the box (inverse of flexibility)
    K = np.linalg.inv(F)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Axial Properties ============================================================
    
    
    # Axial stiffness matrix A
    A = np.zeros((3,3))
    A[ 0 , 0 ] = K[0  , 0 ]
    A[0  ,1:3] = K[0  ,4:6]
    A[1:3, 0 ] = K[4:6, 0 ] 
    A[1:3,1:3] = K[4:6,4:6] 
    
    
    # Axial center
    
    # Axial center coordinates
    y_NA = -A[0,2]/A[0,0]
    z_NA = +A[0,1]/A[0,0]
    
    
    # Translation matrix to uncouple axial force and bending moments
    T_axial = np.array([[  1  , 0, 0], 
                        [-z_NA, 1, 0],
                        [+y_NA, 0, 1]])
    
    # A in the translated SRef
    AA = np.dot( np.dot(T_axial,A) , np.transpose(T_axial) )
    
    
    # Rotation angle that decouples the bending moments
    # In case they are already uncoupled (AA23=0 already) alpha angle should be zero
    if AA[1,2]==0:
        alpha=0.0
    else:
        alpha = 1/2 * np.arctan( 2*AA[1,2] / (AA[1,1]-AA[2,2]) )
    
    
    
    # Rotation matrix
    R_axial = np.array([[     1     ,        0       , 0             ],
                        [     0     ,  np.cos(alpha) , np.sin(alpha) ], 
                        [     0     , -np.sin(alpha) , np.cos(alpha) ]])
    
    # A in the translated and rotated SRef
    AAA = np.dot( np.dot(R_axial,AA) , np.transpose(R_axial) )
    
    
    
    # Properties identification
    EA  =  AA[0,0] #=AAA[0,0]
    EIy = AAA[1,1]
    EIz = AAA[2,2]
    
    
    a  = EA /E
    Iy = EIy/E
    Iz = EIz/E 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Displacement correction =====================================================
    
    # A correction involving displacements uy,uz due to unit forces Fy,Fz is needed
    # These displacements are the result of a component due to shear and another
    # due to bending
    # For the calculation of the shear properties, we want only the contribution
    # due to shear
    # Note: this affectes elements D11 and D22 of the flexibility matrix
    # which have not been used previosly for axial properties
    # Furthermore, Iy,Iz are needed, that s why we apply now the correction
    # Eq.(3.25) Cirilo: u_shear = u - u_bending
    # where u_bending is computed from an analytical formula
    
    # D22 (D11 in python) -> displacement uy
    # Fy produces uy due to bending also, not only due to shear
    # D33 (D22 in python) -> displacement uz
    # Fz produces uz due to bending also, not only due to shear
    D1[1,1] = D1[1,1] - l1**3/(3*EIz) + (l2-l1)*l1**2/(2*EIz)
    D1[1,1] = D1[2,2] - l1**3/(3*EIy) + (l2-l1)*l1**2/(2*EIy)
    D2[1,1] = D2[1,1] - l2**3/(3*EIz)
    D2[2,2] = D2[2,2] - l2**3/(3*EIy)
    
    
    F = (D2-D1)/(l2-l1)
    
    F[ 0  , 1:4] = 0 #row  1,    , columns 2,3,4
    F[1:4 ,  0 ] = 0 #rows 2,3,4 , column  1
    F[1:4 , 4:6] = 0 #rows 2,3,4 , columns 5,6
    F[4:6 , 1:4] = 0 #rows 5,6   , columns 2,3,4
    
    # Modified Stiffness Matrix of the box
    K = np.linalg.inv(F)
    
    
    
    # Correction seems not to affect shear center of J values
    # And Kz Ky are not well predicted anyway
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Shear Properties ============================================================
    
    
    
    # Shear stiffness matrix S
    S = np.zeros((3,3))
    S[0:3,0:3] = K[1:4,1:4]
    
    
    # Rotation angle that decouples the shear forces
    # In case they are already uncoupled (S12=0 already) beta angle should be zero
    if S[0,1]==0:
        beta=0.0
    else:    
        beta = 1/2 * np.arctan( -2*S[0,1] / (S[1,1]-S[0,0]) )
    
    
    
    # Rotation matrix to uncouple the shear forces
    R_shear = np.array([[  np.cos(beta) ,  np.sin(beta), 0], 
                        [ -np.sin(beta) ,  np.cos(beta), 0],
                        [      0        ,     0        , 1]])
               
    
    # S in the rotated SRef
    SS = np.dot( np.dot(R_shear,S) , np.transpose(R_shear) )
    
    
    # Identification of section properties
    GAy = SS[0,0]
    GAz = SS[1,1]   
    
    
    # Shear centre coordinates
    y_SC =  SS[1,2]/GAz 
    z_SC = -SS[0,2]/GAy 
    
    
    # Translation matrix to uncouple shear forces and torsion
    T_shear = np.array([[  1  ,   0   , 0 ], 
                        [  0  ,   1   , 0 ],
                        [ z_SC, -y_SC , 1 ]])
    
    # S in the rotated and translated SRef
    SSS = np.dot( np.dot(T_shear,SS) , np.transpose(T_shear) )
    
    # Section properties
    GJ = SSS[2,2]
    J  = GJ/G
    Ky = abs( GAy/(G*a) )
    Kz = abs( GAz/(G*a) )
    
    
    
    
    
    
    
    
    # Results =====================================================================
    #alpha_deg = alpha*180/np.pi
    #beta_deg  =  beta*180/np.pi
    results = {'y_NA':y_NA, 'z_NA':z_NA, 'alpha':alpha, 'a':a, 'Iy':Iy, 'Iz':Iz,\
    'beta':beta, 'y_SC':y_SC, 'z_SC':z_SC, 'J':J, 'Ky':Ky, 'Kz':Kz}    
    
    return(results)    
    
   












