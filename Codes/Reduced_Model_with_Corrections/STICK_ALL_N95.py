# -*- coding: utf-8 -*-


# Input file:
    # Nastran_template.inp

# Output files:
    # FEM input (.inp) and results (.pnh) files of 3D (full) and reduced (stick) models
    # FullStaticM .inp and .pnh
    # FullModal   .inp and .pnh
    # StickStaticM.inp and .pnh
    # StickModal  .inp and .pnh

#  Post-Processing results:
    # Static response comparison from unit loads at the tip 
    # Natural frequencies and mode shapes comparsion (MAC)



# *****************************************************************************
# *****************************************************************************
# ********************  STICK_STATIC    ***************************************
# *****************************************************************************
# *****************************************************************************

print('\n'*50+" ======== Execution ======== \n")

import numpy as np
import nastran_external_code_v2
import All_Functions
import writing_cards_N95 as cards
from All_Functions import free as free
import Extraction_Functions
import os
import geometry_v2


# =============================================================================
# Inputs (the ones that can change in the optimization process)
# =============================================================================
# SI units
# Geometrical Inputs
ChordFactor = 1
ScaleFactor = 26.2
# FEM inputs
#Vector containing the thickness of each region
# [ribs, skin1-skin4,sparLE1-sparLE4,sparMID1-sparMID4,sparTE1-sparTE4]  1=root,4=tip
t = [0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010,0.010]
#Vector containing the cross section of the stringers
s = [0.000342]
#Young's modulus
E = 70.E+9
#Poisson's ratio
nu = 0.33
#Material density
rho_s = 2700.
# Required file as input: Nastran_template.inp (Nastran template)










# *****************************************************************************
# SECTION 1: GEOMETRY GENERATION
#   Nastran_template.inp
# outputs:
#   full_basic_input.inp
#   struct.igs
#   geometry.geo
#   struc_mesh.bdf
# *****************************************************************************



# =============================================================================
# Get coodinates
# =============================================================================
# We generate the geometry and mesh, we read coordinates to fill the template
geometry_object = geometry_v2.Geometry(ScaleFactor, ChordFactor)       
node_coord_all  = geometry_object.solve_nonlinear()['node_coord_all']
#struct_output_data['node_coord_all'][0][1]

# =============================================================================
# Modify template
# =============================================================================
# We change variables in the template by their values
nastran_object = nastran_external_code_v2.Nastran(node_coord_all, t, s, E, nu, rho_s)
nastran_object.solve_nonlinear()







   





# *****************************************************************************
# SECTION 2: EXTRACTION OF EQUIVALENT BEAM SECTION (STIFFNESS) PROPERTIES
# inputs:
    # struct_mesh.bdf
    # full_basic.inp 
# outputs:
    # section properties (dictionary S_P) 
# *****************************************************************************

  
# =============================================================================
# Full Model Basic Nastran Input
# =============================================================================
# The filled template provides a basic nastran input
# with two markers, where we will add the things we need for each study case
#$ADD_PRE_BULK_CARDS
#$ADD_BULK_CARDS
file_basic_input  = open('full_basic.inp','r')
text_basic_input  = file_basic_input.read()
lines_basic_input = text_basic_input.splitlines()
file_basic_input.close() 


# Some things depend on the template used (marked as *template_specific)


# =============================================================================
# Full Model Unit Load Static Analyses for Properties Extraction Nastran Input
# =============================================================================
full_model_extraction_name = 'FullStatic.inp'
central_id = 7000 # Initial id of Central nodes and RBE2 elements
SRef_id    = 77
n_ribs     = 21 #*template_specific

# Central nodes ids and coordinates ===
# ids of central nodes
centrals = central_id*np.ones((n_ribs,), dtype=np.int) + range(0,n_ribs)
LE_rib_corner_ids    = range(1,201+1,10) #*template_specific
LE_rib_corner_coords = All_Functions.grid_coordinates('struct_mesh.bdf', LE_rib_corner_ids)
ribs_y_coords = []
for node_id in LE_rib_corner_ids:
    ribs_y_coords.append(LE_rib_corner_coords[node_id][1])
# Stick model nodes (called centrals) in stick model frame
centrals_Xcoords = ribs_y_coords

# Stick model origin ===
# Stick model axis is place parallel to geometry's y-axis at the location of the 
# middle spar part which is parallel to this axis (central spar part away from the wing root)
# X Coordinates of stick model frame origin wrt geometrical frame origin
xgeom_stick_origin = All_Functions.grid_coordinates('struct_mesh.bdf', [476])[476][0] #*template_specific

# RBE2 dependent nodes ===
# List of lists of nodes to which we link the central nodes with RBE2 at each section
nodes_sections = []
for nn in range(0,n_ribs):
    nodes_sections.append([nn*10+1,nn*10+4,nn*10+5,nn*10+6,nn*10+9,nn*10+10]) 

# Bdf writing ===
full_model_extraction_file = open(full_model_extraction_name,'w')
stop = 0
for i,line in enumerate(lines_basic_input):
    if '$ADD_PRE_BULK_CARDS' in line:
        full_model_extraction_file.write(line+'\n')
        full_model_extraction_file.write('ID,CRM Wing Statics\n')
        full_model_extraction_file.write('SOL,1\n')
        full_model_extraction_file.write('CEND\n')
        full_model_extraction_file.write('TITLE = CRM_Wing Full Model Static Analyses For Extraction\n')
        full_model_extraction_file.write('SPC=101\n')
        full_model_extraction_file.write('$Marker to place python generated subcases\n')
        full_model_extraction_file.write('$${\n$$}\n')
    elif '$ADD_BULK_CARDS' in line:
        full_model_extraction_file.write(line+'\n')
        full_model_extraction_file.write('$Central nodes\n')
        for i,id in enumerate(centrals):
            full_model_extraction_file.write(cards.GRID(id,SRef_id,[free(centrals_Xcoords[i]), 0., 0.],SRef_id))
        full_model_extraction_file.write('$RBE2 elements (CRGID1)\n')
        for i,node in enumerate(centrals):
            full_model_extraction_file.write(cards.CRIGD1(central_id+i,centrals[i],nodes_sections[i]))
        full_model_extraction_file.write('$Reference System for the analysis (x=beam axis)\n')
        full_model_extraction_file.write(cards.CORD2R(SRef_id,[free(xgeom_stick_origin),0.,0.],[free(xgeom_stick_origin),0.,1.],[free(xgeom_stick_origin),1.,0.]))
        full_model_extraction_file.write('$Marker to place python generated loads\n')
        full_model_extraction_file.write('$${\n$$}\n')
        full_model_extraction_file.write('ENDDATA\n')
    elif 'GRID' in line:
        line += ','+str(SRef_id)
        full_model_extraction_file.write(line+'\n')
    # Modification of the boundary conditions, beacause of the RBE2
    elif 'SPC1' in line:
        line+=str(centrals[0])+','
        spc_lines = [line,lines_basic_input[i+1],lines_basic_input[i+2],lines_basic_input[i+3]]
        for spc_line in spc_lines:
            for item in nodes_sections[0]:
                target = ','+str(item)+','
                if target in spc_line:
                    spc_line = spc_line.replace(target,',')
            full_model_extraction_file.write(spc_line+'\n')
        stop = 1
    elif not stop:
        full_model_extraction_file.write(line+'\n')
full_model_extraction_file.close()

# UNIT LOAD ANALYSES FOR THE SECTION PROPERTIES EXTRACTION ====================
full_model_extraction_name_mod = Extraction_Functions.generate_loads(list(centrals),full_model_extraction_name,SRef_id)




# =============================================================================
# Launch Nastran for static unit load cases for extraction (punch file is generated)
# =============================================================================
# RUN ANALYSES ON WINDOWS AND COPY PASTE PUNCH FILES
# Run Nastran, all unit load analyses are run
# Punch with all displacements is created
print('nastran full_model_extraction_name_mod')
os.system('nastran '+full_model_extraction_name_mod[0:-4])
name_punch = full_model_extraction_name_mod[:-3]+'pnh'





                               

# =============================================================================
# EQUIVALENT BEAM SECTION PROPERTIES EXTRACTION (STIFFNESS)
# =============================================================================
# Inputs ===
G = E/(2*(1+nu))
# coord_grid is also used in the stick model bdf writing
coord_grid = All_Functions.grid_coordinates(full_model_extraction_name, centrals)
MID        = 1 # Material ID

# Computations ===
# Number of boxes
# In normal notation (start at 1) box n goes from central n to central n+1
# e.g. box1: central1-central2
n_centrals = len(centrals)
n_boxes    = n_centrals - 1
boxes      = list(range(1,n_boxes+1))
# Dictionary in which results will be stored
S_P = {}
# Extraction process for each box
for n_box in boxes:
    # Nodes associated to the box
    node_i = centrals[n_box-1]
    node_j = centrals[ n_box ]
    # Subcases associated to the box
    # Each box contains 6 associated subcases (6 unit loads at node_j)
    n_ini = 1 + (n_box-1)*6
    n_fin = n_ini + 6
    subcases = []
    for n_subcase in range(n_ini,n_fin):
        subcases.append('SUBCASE'+str(n_subcase))
    # Flexibility matrices relating loads at nodej and displacements at
    # node_i (D1) and at node_j D2. Obtained from the displacements due to unit loads
    [Di,Dj] = Extraction_Functions.fleximatrices(name_punch, node_i, node_j, subcases)
    # Distance to the clamped end (which is at x=0) (measured in x-direction as x=beam axis)
    # i.e. x-coordinates of nodes i and j
    li = coord_grid[node_i][0]
    lj = coord_grid[node_j][0]
    # Equivalent section properties extraction (S_P stands for Section_Properties)
    S_P[n_box] = Extraction_Functions.properties_extraction(Di,Dj,li,lj,E,G)
#print(S_P)










# *****************************************************************************
# SECTION 3: Stick model construction for static analysis
# inputs:
    # stick_basic.inp
    # section properties
# outputs:
    # StickStaticM.inp 
# *****************************************************************************

   
  
# =============================================================================
# STICK MODEL CONSTRUCTION (Basic Nastran Input)
# =============================================================================
stick_model_name = 'stick_basic.inp'
# Boundary Condition: 1st node clamped
[BC1,BC2]=cards.BC_clamped(centrals[0])
datfile=open(stick_model_name,'w')
datfile.write('$ADD_PRE_BULK_CARDS\n')
datfile.write('BEGIN BULK\n')
# No density, mass will be modeled through concetrated masses
datfile.write('\n\n$Material\n')
datfile.write(cards.MAT1(MID,free(E),free(G),0.))
datfile.write('\n$Nodes\n')
for ID in centrals: #we put ID in centrals instead of in coord_grid to have the correct order in the bdf
    X  = coord_grid[ID]
    XX =[free(X[0]), free(X[1]), free(X[2]) ]
    grid_card  = cards.GRID(ID,0,XX,0)
    datfile.write(grid_card)
datfile.write('\n$Beam elements and properties\n')
# Each box will turn into a beam element
for i in boxes:
    EID = i
    PID = i
    GA  = centrals[i-1]
    GB  = centrals[i]
    # Element y-axis direction
    X   = [0.0, free(np.cos(S_P[1]['beta'])), free(np.sin(S_P[1]['beta'])) ]  #sure? 
    # Shear center position    
    W   = [0.0, free(S_P[i]['y_SC']), free(S_P[i]['z_SC']) ]                  
    # Section Properties    
    A   = free(S_P[i]['a'])
    Iy  = free(S_P[i]['Iy']) 
    Iz  = free(S_P[i]['Iz'])
    J   = free(S_P[i]['J'])
    # Neutral axis position wrt shear center
    [NA_1,NA_2] = [ S_P[i]['y_NA']-S_P[i]['y_SC'],   S_P[i]['z_NA']-S_P[i]['z_SC']  ]
    NA  = ["%.2f"%NA_1,  "%.2f"%NA_2]       
    # CBEAM cards =================================================================
    # EID   -> element id
    # PID   -> property id
    # GA,GB -> the two grid points of the beam element
    # X     -> orientation vector of local y-axis (x=local beam axis)
    # W     -> offset vector of the shear center
    cbar_card = cards.CBAR(EID,PID,GA,GB,X,W)
    # PBEAM cards =================================================================
    # PID           -> propery id
    # MID           -> material id
    # A,I1,I2,I12,J -> section properties (Iyz=0)
    # K1,K2 -> shear stiffness factor are null, i.e. Bernoulli-Euler beam
    # N -> neutral axis coordinates from shear center
    # Same properties all along the beam element
    pbar_card = cards.PBAR(PID,MID,A,Iy,Iz,J) 
    # Writing on the bdf ==========================================================
    #datfile.write('\n$= Beam Element\n')
    datfile.write(cbar_card+'\n')
    #datfile.write('\n$= Beam Property\n')
    datfile.write(pbar_card+'\n')


datfile.write('\n$Boundary conditions\n')
datfile.write(BC2)
datfile.write('$ADD_BULK_CARDS\n')
datfile.write('\nENDDATA\n')
datfile.close()

# We read the stick model basic input
file_stick_basic  = open(stick_model_name,'r')
text_stick_basic  = file_stick_basic.read()
lines_stick_basic = text_stick_basic.splitlines()
file_stick_basic.close() 








  


# =============================================================================
# STICK MODEL unit loads at tip for static comparison
# =============================================================================
stick_model_static = open('StickStatic.inp','w')
for line in lines_stick_basic:
    if '$ADD_PRE_BULK_CARDS' in line:
        stick_model_static.write(line+'\n')
        stick_model_static.write('ID,Stick Statics\n')
        stick_model_static.write('SOL,1\n')
        stick_model_static.write('CEND\n')
        stick_model_static.write('TITLE = Stick Model Static Analyses for Comparison\n')
        stick_model_static.write(BC1+'\n')
        stick_model_static.write('$Marker to place python generated subcases\n')
        stick_model_static.write('$${\n$$}\n')
    elif '$ADD_BULK_CARDS' in line:
        stick_model_static.write('$Marker to place python generated loads\n')
        stick_model_static.write('$${\n$$}\n')
    else:
        stick_model_static.write(line+'\n')
stick_model_static.close()

# UNIT LOAD ANALYSES FOR THE SECTION PROPERTIES EXTRACTION ====================
stick_model_static_name_mod = Extraction_Functions.generate_loads(list(centrals),'StickStatic.inp',0)



 
# =============================================================================
# STICK MODEL unit loads at tip for static comparison
# =============================================================================  
print('nastran stick_model_static_name_mod')  
os.system('nastran '+stick_model_static_name_mod[0:-4])





























# *****************************************************************************
# *****************************************************************************
# ******************** STICK_MODAL    ****************************************
# *****************************************************************************
# *****************************************************************************

import Mass_Functions


# In terms of files:
# Inputs:
# - FullStatic.inp
# - stick_basic.inp
# - full_basic.inp
# - struct_mesh.bdf

# Outputs:
# StickModal .inp and .pnh
# FullModal .inp and .pnh





# *****************************************************************************
# SECTION 4: Extraction of Mass Properties
# inputs:
    # FullStatic.inp
# outputs:
    # Mass Properties (dictionary M_P)
# *****************************************************************************





# =============================================================================
# Ids for Box Creations   *template_specific
# =============================================================================
n_ribs      = 21
n_boxes     = n_ribs-1
n_ini       = 4049
# ids of spar caps and stringers
elem_ids = {}
for j in range(0,5):#3 sets of spar caps, 2 of stringers
    elem_ids['lines'+str(j+1)]=[]
    for i in range(0,n_boxes):
        n_end = n_ini+4
        elem_ids['lines'+str(j+1)].append(range(n_ini,n_end))
        n_ini = n_end       
# ids of the ribs
for j in range(0,n_ribs):
    n_end = n_ini+16
    elem_ids['rib'+str(j+1)]=range(n_ini,n_end)
    n_ini = n_end
# ids of the skin
for j in range(0,n_boxes):
    n_end = n_ini+32
    elem_ids['skin'+str(j+1)]=range(n_ini,n_end)
    n_ini = n_end   
# ids of spars
for j in range(0,3):#3 spars
    elem_ids['spar'+str(j+1)]=[]
    for i in range(0,n_boxes):
        n_end = n_ini+4
        elem_ids['spar'+str(j+1)].append(range(n_ini,n_end))
        n_ini = n_end    
box_elems = {}    
for i in range(0,n_boxes):
    if i==0:
        box_elems['box'+str(i+1)] = elem_ids['lines1'][i] + elem_ids['lines2'][i] + elem_ids['lines3'][i] \
                    + elem_ids['lines4'][i] + elem_ids['lines5'][i] + elem_ids['rib'+str(i+2)] \
                    + elem_ids['skin'+str(i+1)] + elem_ids['spar1'][i] + elem_ids['spar2'][i] + elem_ids['spar3'][i] \
                    + elem_ids['rib'+str(i+1)]
    else:
        box_elems['box'+str(i+1)] = elem_ids['lines1'][i] + elem_ids['lines2'][i] + elem_ids['lines3'][i] \
                    + elem_ids['lines4'][i] + elem_ids['lines5'][i] + elem_ids['rib'+str(i+2)] \
                    + elem_ids['skin'+str(i+1)] + elem_ids['spar1'][i] + elem_ids['spar2'][i] + elem_ids['spar3'][i]             
grid_ids = {}
grid_ids['rib_1_1']          = np.array([1,4,2,3,5,6,7,8,9,10])
grid_ids['rib_1_2']          = np.array([211,213,215,217,218,220,221,223,214,212,216,219,222])
grid_ids['rib_1_3']          = np.array([684,685,686,687])
grid_ids['between_ribs_1_1'] = np.array([486,484,487,485,489,488,491,490,493,492])
grid_ids['between_ribs_1_2'] = np.array([769,768,771,770,773,772,775,774])    
grid_ids['between_ribs_1_3'] = np.array([928,948,968])           
for i in range(0,n_ribs):
    grid_ids['rib_'+str(i+2)+'_1']          = grid_ids['rib_'+str(i+1)+'_1'] + 10
    grid_ids['rib_'+str(i+2)+'_2']          = grid_ids['rib_'+str(i+1)+'_2'] + 13
    grid_ids['rib_'+str(i+2)+'_3']          = grid_ids['rib_'+str(i+1)+'_3'] + 4
    if i!=n_ribs-1:
        grid_ids['between_ribs_'+str(i+2)+'_1'] = grid_ids['between_ribs_'+str(i+1)+'_1'] + 10
        grid_ids['between_ribs_'+str(i+2)+'_2'] = grid_ids['between_ribs_'+str(i+1)+'_2'] + 8
        grid_ids['between_ribs_'+str(i+2)+'_3'] = grid_ids['between_ribs_'+str(i+1)+'_3'] + 1             
for i in range(0,n_ribs):
    grid_ids['rib_'+str(i+1)] = list(grid_ids['rib_'+str(i+1)+'_1']) + list(grid_ids['rib_'+str(i+1)+'_2']) + list(grid_ids['rib_'+str(i+1)+'_3'])
    if i!=n_ribs-1:
        grid_ids['between_ribs_'+str(i+1)] = list(grid_ids['between_ribs_'+str(i+1)+'_1']) + \
                 list(grid_ids['between_ribs_'+str(i+1)+'_2']) + list(grid_ids['between_ribs_'+str(i+1)+'_3'])
            
             







                    
                    
                    
                             
          
# =============================================================================
# CODE CREATING BOX MODELS 
# =============================================================================
# Inputs
# section_gridpoints: list of lists with node ids at each section
# section_elements:   list of lists with elements in each box

# Function creating the box model, if the nodes and elements' ids are given
name_full_model  =  'FullStatic.inp'
boxes = range(1,n_ribs)
box_models  = []
box_results = []
# A model is created for each box
for box in boxes:
    the_box_grids  = grid_ids['rib_'+str(box)] + grid_ids['between_ribs_'+str(box)] + grid_ids['rib_'+str(box+1)]
    the_box_elems  = box_elems['box'+str(box)]
    name_box_model   = 'Box_Model_'+str(box)+'.inp'
    name_box_result  = 'Box_Model_'+str(box)+'.out'
    box_models.append(name_box_model)
    box_results.append(name_box_result)
    Mass_Functions.box_model(the_box_grids,the_box_elems,name_full_model,name_box_model)      


# --- Bookmark 1     This will appear in the outline explorer 
# =============================================================================
# CODE LAUNCHING NASTRAN ANALYSIS (JUST TO GET MASS PROPERTIES)
# =============================================================================   
#==============================================================================
# for box in boxes:
#    os.system('nastran '+'Box_Model_'+str(box))
#    print('Nastran is executing Box Model '+str(box))
# 
#==============================================================================


# =============================================================================
# CODE CREATING THE MASS PROPERTIES
# ============================================================================= 
# Read mass properties from the .out of each box 
M_P = {}
for i,result_file in enumerate(box_results):
    #print(result_file)
    [MASS,XCG,YCG,ZCG,IXX,IYY,IZZ,IXY,IYZ,IXZ] = Mass_Functions.read_mass_properties(result_file)
    M_P['Box'+str(i+1)] = {}
    M_P['Box'+str(i+1)]['mass'] = MASS
    M_P['Box'+str(i+1)]['xcg']  = XCG - ScaleFactor/n_boxes*i
    M_P['Box'+str(i+1)]['ycg']  = YCG
    M_P['Box'+str(i+1)]['zcg']  = ZCG
    M_P['Box'+str(i+1)]['Ixx']  = IXX
    M_P['Box'+str(i+1)]['Iyy']  = IYY
    M_P['Box'+str(i+1)]['Izz']  = IZZ
    M_P['Box'+str(i+1)]['Ixy']  = IXY
    M_P['Box'+str(i+1)]['Iyz']  = IYZ
    M_P['Box'+str(i+1)]['Ixz']  = IXZ
    
#print(M_P)
#print('\n')
#print(M_P['Box1']['ycg'])







# *****************************************************************************
# SECTION 5: Creation and execution of stick and full models for modal analysis
# inputs:
    # struct_mesh.bdf
    # Mass Properties
    # full_basic.inp
    #stick_basic.inp
# outputs:
    # FullModal .inp and .pnh
    # StickModal .inp and .pnh
# *****************************************************************************



# =============================================================================
# STICK MODEL for Modal Analysis (with mass properties)
# =============================================================================
central_id = 7000 # Initial id of Central nodes and RBE2 elements (same that in STICK_STATIC.py code)
centrals = central_id*np.ones((n_ribs,), dtype=np.int) + range(0,n_ribs)
mass_nodes = centrals[0:-1]
SRef_id = 77
# ids of nodes used for mode shape comparison
stick_LE_ids = range(214,474+1,13)#*template_specific
stick_TE_ids = range(222,482+1,13)#*template_specific  
# reading coordinates  for mode shape comparison nodes
stick_LE_coords = All_Functions.grid_coordinates('struct_mesh.bdf', stick_LE_ids)
stick_TE_coords = All_Functions.grid_coordinates('struct_mesh.bdf', stick_TE_ids)

# Stick model origin ===
# Stick model axis is placed parallel to geometry's y-axis at the location of the 
# middle spar part which is parallel to this axis (central spar part away from the wing root)
# X Coordinates of stick model frame origin wrt geometrical frame origin
xgeom_stick_origin = All_Functions.grid_coordinates('struct_mesh.bdf', [476])[476][0] #*template_specific 


# We read the stick model basic input
file_stick_basic  = open('stick_basic.inp','r')
text_stick_basic  = file_stick_basic.read()
lines_stick_basic = text_stick_basic.splitlines()
file_stick_basic.close()

stick_model_modal_name = 'StickModal.inp'
stick_model_modal = open(stick_model_modal_name,'w')
for line in lines_stick_basic:
    if '$ADD_PRE_BULK_CARDS' in line:
        stick_model_modal.write(line+'\n')
        stick_model_modal.write('ID STICK MODAL\n')
        stick_model_modal.write('SOL,3\n')
        stick_model_modal.write('CEND\n')
        stick_model_modal.write('TITLE = Stick Model Modal Analysis n\n')
        stick_model_modal.write('METHOD = 1\n')
        stick_model_modal.write('SPC = 1\n')
        stick_model_modal.write('SET 1 = ')
        for ii,item in enumerate(stick_LE_ids+stick_TE_ids):
            if ii==len(stick_LE_ids+stick_TE_ids)-1:
                stick_model_modal.write(str(item)+'\n')
            elif ii%9==0:
                stick_model_modal.write(str(item)+',\n')
            else:
                stick_model_modal.write(str(item)+',')
        
        stick_model_modal.write('\n\nVECTOR(PUNCH)=1\n')

    elif '$ADD_BULK_CARDS' in line:
        stick_model_modal.write('PARAM,GRDPNT,0\n')
        stick_model_modal.write('EIGR,1,FEER,0.1,,,10,,,\n,MAX\n')
        # Mass properties
        stick_model_modal.write('\n$Concentrated masses\n') 
        
        for node_id in stick_LE_ids:
            stick_LE_x_coord = stick_LE_coords[node_id][1]
            stick_LE_y_coord = float(stick_LE_coords[node_id][0]) - float(xgeom_stick_origin)
            stick_model_modal.write(cards.GRID(node_id,'0',[free(stick_LE_x_coord), free(stick_LE_y_coord), 0.],'0'))
        for node_id in stick_TE_ids:        
            stick_TE_x_coord = stick_TE_coords[node_id][1]
            stick_TE_y_coord = float(stick_TE_coords[node_id][0]) - float(xgeom_stick_origin)
            stick_model_modal.write(cards.GRID(node_id,'0',[free(stick_TE_x_coord), free(stick_TE_y_coord), 0.],'0'))
        for i,id in enumerate(centrals):
            stick_model_modal.write(cards.CRIGD1(8000+i,centrals[i],[stick_LE_ids[i],stick_TE_ids[i]]))
            
            
        for i,node in enumerate(mass_nodes):
            
            EID   = n_boxes+i+1
            G     = node
            BMP   = M_P['Box'+str(boxes[i])] # Box Mass Properties
            M_i   = BMP['mass']/3
            M_i   =free(M_i)
            XCG_i = [free(BMP['xcg']), free(BMP['ycg']), free(BMP['zcg'])] 
            Ixx_i = free(BMP['Ixx']) 
            Iyy_i = free(BMP['Iyy']) 
            Izz_i = free(BMP['Izz']) 
            Ixy_i = free(BMP['Ixy']) 
            Ixz_i = free(BMP['Ixz']) 
            Iyz_i = free(BMP['Iyz'])
            stick_model_modal.write(cards.CONM2(EID,G,M_i,XCG_i,Ixx_i,Iyy_i,Izz_i,Ixy_i,Ixz_i,Iyz_i))
        p=40    
        for i in range(20):
            p+=1
            EID   = p
            G     = stick_LE_ids[i]
            BMP   = M_P['Box'+str(boxes[i])] # Box Mass Properties
            M_i   = BMP['mass']/3
            M_i   =free(M_i)
            node_id=G
            
            XCG_i = [0, 0, 0] 
            Ixx_i = free(BMP['Ixx']) 
            Iyy_i = free(BMP['Iyy']) 
            Izz_i = free(BMP['Izz']) 
            Ixy_i = free(BMP['Ixy']) 
            Ixz_i = free(BMP['Ixz']) 
            Iyz_i = free(BMP['Iyz'])
            stick_model_modal.write(cards.CONM2(EID,G,M_i,XCG_i,Ixx_i,Iyy_i,Izz_i,Ixy_i,Ixz_i,Iyz_i))

        for i in range(20):    
            p+=1
            EID   = p
            G     = stick_TE_ids[i]
            BMP   = M_P['Box'+str(boxes[i])] # Box Mass Properties
            M_i   = BMP['mass']/3
            M_i   =free(M_i)
            node_id=G
            
            XCG_i = [0, 0, 0] 
            Ixx_i = free(BMP['Ixx']) 
            Iyy_i = free(BMP['Iyy'])
            Izz_i = free(BMP['Izz']) 
            Ixy_i = free(BMP['Ixy']) 
            Ixz_i = free(BMP['Ixz']) 
            Iyz_i = free(BMP['Iyz'])
            stick_model_modal.write(cards.CONM2(EID,G,M_i,XCG_i,Ixx_i,Iyy_i,Izz_i,Ixy_i,Ixz_i,Iyz_i))
        stick_model_modal.write('\n$RBE2 for mode shape comparison\n')
        # Reads in 3D global/geometry frame. Writes in stick local frame
        
      
            
    else:
        stick_model_modal.write(line+'\n')
stick_model_modal.close()



# =============================================================================
# Full Model Modal Analysis Nastran Input
# ============================================================================= 

# Reading of the basic nastran input .bdf
file_basic_input  = open('full_basic.inp','r')
text_basic_input  = file_basic_input.read()
lines_basic_input = text_basic_input.splitlines()
file_basic_input.close() 
xgeom_stick_origin = All_Functions.grid_coordinates('struct_mesh.bdf', [476])[476][0] #*template_specific
# Writing of the full model modal .bdf 
full_model_modal_name = 'FullModal.inp'
full_model_modal_file = open(full_model_modal_name,'w')
for line in lines_basic_input:
    if '$ADD_PRE_BULK_CARDS' in line:
        full_model_modal_file.write(line+'\n')
        full_model_modal_file.write('ID STICK MODAL\n')
        full_model_modal_file.write('SOL,3\n')
        full_model_modal_file.write('CEND\n')
        full_model_modal_file.write('TITLE = CRM_Wing Full Model Modal\n')
        full_model_modal_file.write('METHOD=1\n')
        full_model_modal_file.write('SPC=101\n')          
        full_model_modal_file.write('SET 1 = ')
        for ii,item in enumerate(stick_LE_ids+stick_TE_ids):
            if ii==len(stick_LE_ids+stick_TE_ids)-1:
                full_model_modal_file.write(str(item)+'\n')
            elif ii%9==0:
                full_model_modal_file.write(str(item)+',\n')
            else:
                full_model_modal_file.write(str(item)+',')   
        full_model_modal_file.write('\n\nVECTOR(PUNCH)=1\n')
    elif '$ADD_BULK_CARDS' in line:
        full_model_modal_file.write(line+'\n')
        full_model_modal_file.write('PARAM,GRDPNT,0\n')
        full_model_modal_file.write('EIGR,1,FEER,0.1,,,10,,,\n,MAX\n')
        full_model_modal_file.write('$Reference System for the analysis (x=beam axis)\n')
        full_model_modal_file.write(cards.CORD2R(SRef_id,[free(xgeom_stick_origin),0.,0.],[free(xgeom_stick_origin),0.,1.],[free(xgeom_stick_origin),1.,0.]))
    elif 'GRID' in line:
        line += ','+str(SRef_id)
        full_model_modal_file.write(line+'\n')
    else:
        full_model_modal_file.write(line+'\n')
full_model_modal_file.close()

# =============================================================================
# Launch modal analysis of both full and stick models
# ============================================================================= 
print('nastran modal full')
os.system('nastran '+full_model_modal_name[0:-4])  # full  Model
print('nastran modal stick')
os.system('nastran '+stick_model_modal_name[0:-4]) # stick Model










         
         
         
         
# *****************************************************************************
# *****************************************************************************
# ********************  POST_PROCESSING STATICS    ****************************
# *****************************************************************************
# *****************************************************************************         
             
from Extraction_Functions import read_punch
import matplotlib.pyplot as plt

# Punch files with displacements due to unit loads at the tip
punch1 = 'FullStaticM.pnh'
punch2 = 'StickStaticM.pnh'

# Read the displacements from the punch file
Disp1 = read_punch(punch1)
Disp2 = read_punch(punch2)

# Subcases corresponding to tip loads
subcases1 = ['SUBCASE115','SUBCASE116','SUBCASE117',\
'SUBCASE118','SUBCASE119','SUBCASE120']
subcases2 = ['SUBCASE115','SUBCASE116','SUBCASE117',\
'SUBCASE118','SUBCASE119','SUBCASE120']

#central_id = 7000 # Initial id of Central nodes and RBE2 elements
#n_ribs     = 21
## ids of central nodes
#centrals = central_id*np.ones((n_ribs,), dtype=np.int) + range(0,n_ribs)

# we need the nodes where displacements are measured
nodes1  = list(centrals)
Nnodes1 = len(nodes1)
nodes2  = list(centrals)
Nnodes2 = len(nodes2)

## we need the x-coordinates
x_coords = centrals_Xcoords #To have this, run STICK_STATIC.py first

# We will put the displacements in vectors
TX1 = np.zeros(Nnodes1)
TY1 = np.zeros(Nnodes1)
TZ1 = np.zeros(Nnodes1)
RX1 = np.zeros(Nnodes1)
RY1 = np.zeros(Nnodes1)
RZ1 = np.zeros(Nnodes1)
TX2 = np.zeros(Nnodes2)
TY2 = np.zeros(Nnodes2)
TZ2 = np.zeros(Nnodes2)
RX2 = np.zeros(Nnodes2)
RY2 = np.zeros(Nnodes2)
RZ2 = np.zeros(Nnodes2)


# We extract only primary displacements,
# x-disp due to x-force, y-disp due to y-force, etc
for i,node in enumerate(nodes1):
    TX1[i]=Disp1[subcases1[0]][node]['TX']
    TY1[i]=Disp1[subcases1[1]][node]['TY']
    TZ1[i]=Disp1[subcases1[2]][node]['TZ']
    RX1[i]=Disp1[subcases1[3]][node]['RX']
    RY1[i]=Disp1[subcases1[4]][node]['RY']
    RZ1[i]=Disp1[subcases1[5]][node]['RZ']

for i,node in enumerate(nodes2):
    TX2[i]=Disp2[subcases2[0]][node]['TX']
    TY2[i]=Disp2[subcases2[1]][node]['TY']
    TZ2[i]=Disp2[subcases2[2]][node]['TZ']
    RX2[i]=Disp2[subcases2[3]][node]['RX']
    RY2[i]=Disp2[subcases2[4]][node]['RY']
    RZ2[i]=Disp2[subcases2[5]][node]['RZ']



# Ploting ==========================================
#==============================================================================
# LL = ScaleFactor
# x_L = []#non dimensional x coordinate x/L
# for item in x_coords:
#     x_L.append(item/LL)
# 
# def ploting(x,y1,y2,ylab,load): 
#     plt.plot(x, y1, 'bo-',label="full")
#     plt.plot(x, y2, 'gv-',label="stick_cbeam")
#     plt.xlabel('$x/L$', fontsize=14,)
#     plt.ylabel(ylab+' due to unitary '+load, fontsize=14)
#     plt.legend()
#     plt.grid()
#     plt.xlim([0,1])
#     plt.show()
#     
# ploting(x_L,TX1,TX2,'$x$-displ. [m]','$F_x$')
# ploting(x_L,TY1,TY2,'$y$-displ. [m]','$F_y$')
# ploting(x_L,TZ1,TZ2,'$z$-displ. [m]','$F_z$')
# ploting(x_L,RX1,RX2,'$x$-rot. [rad]','$M_x$')
# ploting(x_L,RY1,RY2,'$y$-rot. [rad]','$M_y$')
# ploting(x_L,RZ1,RZ2,'$z$-rot. [rad]','$M_z$')
# 
# 
# 
# 
#    
#==============================================================================
         



# *****************************************************************************
# *****************************************************************************
# ********************  POST_PROCESSING MODAL    ******************************
# *****************************************************************************
# *****************************************************************************  

import Modal_Comparison_Functions as funcs

# MODE SHAPES READING =========================================================
eigs_stick, mode_shapes_stick = funcs.read_punch_modal('StickModal.pnh')
eigs_full , mode_shapes_full  = funcs.read_punch_modal('FullModal.pnh')

N_modes = 10
modes = []
for ii in range(1,N_modes+1):
    modes.append('MODE'+str(ii))
#dofs  = ['TX','TY','TZ','RX','RY','RZ']
dofs  = ['TX','TY','TZ']
#dofs  = ['RX','RY','RZ']
mode_shapes_stick_stacked = {}
mode_shapes_full_stacked  = {}
for mode in modes:
    mode_shapes_stick_stacked[mode] = []
    mode_shapes_full_stacked[mode]  = []
    for dof in dofs:
        mode_shapes_stick_stacked[mode] += mode_shapes_stick[mode][dof]
        mode_shapes_full_stacked[mode]  += mode_shapes_full[mode][dof]
        
MACs=np.zeros((N_modes,N_modes))
funcs.MAC_plot(mode_shapes_stick_stacked,mode_shapes_full_stacked,modes,N_modes,MACs)

# MODE MATCHING: TO BE DONE MANUALLY """"" NOT NOW""""
# From observation of MAC matrix and maybe also of the mode shape in patran
# we match modes:
# Full    Stick
#   1  =>  1
#   2  =>  2
#   3  =>  3
#   4  =>  6
#   6  =>  5
#==============================================================================
# ORDER MODES FULL-STICK
#==============================================================================

match_full  = np.arange(1, N_modes+1)
match_stick = np.zeros(N_modes,dtype=np.int)
for i in range (N_modes):
    for j in range (N_modes):
        if MACs[j,i]>0.75:
            match_stick[i]=j+1

#==============================================================================
# for t in range (N_modes):
#     if match_stick[t]==0:
#         dofs  = ['RX','RY','RZ']
#         
# 
# for mode in modes:
#     mode_shapes_stick_stacked[mode] = []
#     mode_shapes_full_stacked[mode]  = []
#     for dof in dofs:
#         mode_shapes_stick_stacked[mode] += mode_shapes_stick[mode][dof]
#         mode_shapes_full_stacked[mode]  += mode_shapes_full[mode][dof]
# print(match_full)  
# print(match_stick)
# MACs=np.zeros((N_modes,N_modes))
# funcs.MAC_plot(mode_shapes_stick_stacked,mode_shapes_full_stacked,modes,N_modes,MACs)
# for t in range (N_modes):
#     if match_stick[t]==0:
#          for j in range (N_modes):
#              if MACs[j,t]>0.75:
#                  match_stick[t]=j+1
# print(match_full)  
# print(match_stick)
#==============================================================================
remove=np.array([]) 
for i in range (N_modes):
    if match_stick[i]==0:
        remove=np.append(remove,i)
for i in range(len(remove)):
    
    match_full=np.delete(match_full,remove[len(remove)-i-1])
    match_stick=np.delete(match_stick,remove[len(remove)-i-1])
       
print(match_full)  
print(match_stick) 

# Frequencies in Hz
eigs_stick  = map(float,eigs_stick)
eigs_full   = map(float,eigs_full)
freqs_stick = np.sqrt(eigs_stick)/(2*np.pi)
freqs_full  = np.sqrt(eigs_full)/(2*np.pi)     
# Frequency error bar plot
f_error = []
bar_x = []
for ii in range(len(match_full)):
    f_stick = freqs_stick[match_stick[ii]-1]
    f_full  = freqs_full[match_full[ii]  -1]
    f_error.append( (f_full-f_stick)/f_full*100 )
    bar_x.append( str(match_full[ii])+' / '+str(match_stick[ii]) )

y_pos = np.arange(len(f_error)) 
plt.bar(y_pos, f_error, align='center', alpha=0.5)
plt.xticks(y_pos, bar_x)
plt.xlabel('Full / Stick  model modes',size=16)
plt.ylabel('Freq. error $(full-stick)/full$ (%)',size=16)
plt.grid()
plt.show()
print f_error



