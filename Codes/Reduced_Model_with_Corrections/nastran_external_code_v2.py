# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 15:56:23 2016

@author: jmascolo
"""

from __future__ import print_function

#from openmdao.api import ExternalCode

import numpy as np

import os.path, time

#from field_writer_8 import print_float_8
from All_Functions import free as print_float_8

#from is_number import isfloat, isint

class Nastran(object):
    template_file = 'Nastran_template.inp'
    
    #output_file = 'nastran_input_modal.out'
    
    def __init__(self, node_coord_all, t, s, E, nu, rho_s):
        super(Nastran, self).__init__()
    
    
        self.input_filepath = 'full_basic.inp'
        
        #self.output_filepath = 'nastran_input.pnh'
    
#    def __init__(self, node_id, node_id_all, n_stress, tn, sn):
#        super(Nastran, self).__init__()
        
#        #Identification number of the outer surface nodes        
#        self.node_id = node_id
        
#        #Number of nodes on the outer surface        
#        self.ns = len(node_id)
        
#        #Number of stress outputs
#        self.n_stress = n_stress
        
        #Total number of structure nodes
        self.ns_all = 987

        #ID of all nodes
        self.node_id_all = map(str,range(1,self.ns_all+1))
        
        #Number of regions where the thicknesses are defined
        self.tn = 17
        
        #Number of regions where the stringers sections are defined
        self.sn = 1
        
        #Dictionary containing the last used Nastran input parameters
        self.old_params = {}
        
#        #Forces on the nodes of the outer surface        
#        self.add_param('f_node', val=np.zeros((self.ns, 3)))
        
        #Coordinates of all structural nodes
        #add_param('node_coord_all', val=np.zeros((self.ns_all, 3)))
        self.node_coord_all = node_coord_all
        
        #Vector containing the thickness of each region
        #self.add_param('t', val=np.zeros(self.tn))
        self.t = t
        
        #Vector containing the cross section of the stringers
        #self.add_param('s', val=np.zeros(self.sn))
        self.s = s
        
        #Young's modulus
        #self.add_param('E', val=1.)
        self.E = E
        
        #Poisson's ratio
        #self.add_param('nu', val=0.3)
        self.nu = nu
        
        #Material density
        #self.add_param('rho_s', val=1.)
        self.rho_s = rho_s
        
#        #Displacements of the nodes on the outer surface        
#        self.add_output('u', val=np.zeros((self.ns, 3)))
        
#        #Von Mises stress of all elements
#        self.add_output('VMStress', val=np.zeros(self.n_stress))
        
#        #Structural mass
#        self.add_output('M', val=1.)
        
        
#        #Check if the files exist (optional)        
#        self.options['external_input_files'] = [self.input_filepath,]
#        #self.options['external_output_files'] = [self.output_filepath,]
#        
#        self.options['command'] = ['cmd.exe', '/c', r'nastran.bat', self.input_filepath.rstrip('.inp')]
        
#    def solve_nonlinear(self, params, unknowns, resids):   
    def solve_nonlinear(self):
        
        #Copy the Nastran input parameters that affect the stiffnes matrix
#        new_params = {}
#        new_params['node_coord_all'] = params['node_coord_all']
#        new_params['t'] = params['t']
#        new_params['E'] = params['E']
#        new_params['nu'] = params['nu']
#        new_params['rho_s'] = params['rho_s']
                
#        #Compare old and new Nastran input parameters 
#        if new_params == self.old_params:
#            changed = False
#        else:
#            changed = True
        
        # Generate the input file for Nastran from the input file template and pressure values at the nodes
        #The variable "changed" indicates whether the Nastran input parameters have been changed in order to compute or reuse the stiffness matrix
        self.create_input_file()
        
#        # Parent solve_nonlinear function actually runs the external code
#        super(Nastran, self).solve_nonlinear(params, unknowns, resids)
#        
#        output_data = self.get_output_data()
#        
#        # Parse the output file from the external code and set the value of u
#        unknowns['u'] = output_data['u']
#        
#        #Parse the output file from the external code and get the Von Mises Stresses
#        unknowns['VMStress'] = output_data['VMStress']
#        
#        #Parse the output file from the external code and get the structural mass        
#        unknowns['M'] = output_data['M']
#        
#        #Update the old parameters dictionary
#        self.old_params = new_params.copy
  
#    def create_input_file(self, params, changed):       
    def create_input_file(self): 
  
        
#        f_node = params['f_node']
#        node_coord_all = params['node_coord_all']
#        t = params['t']
#        s = params['s']
#        E = params['E']
#        nu = params['nu']
#        rho_s = params['rho_s']
        
        node_coord_all = self.node_coord_all
        t = self.t
        s = self.s
        E = self.E
        nu = self.nu
        rho_s = self.rho_s
        
        input_data = {}
            
#        #Assign each force value to its corresponding node ID in the input data dictionary            
#        for i in range(len(f_node)):
#            input_data['Fx'+self.node_id[i]] = print_float_8(f_node[i, 0])
#            input_data['Fy'+self.node_id[i]] = print_float_8(f_node[i, 1])
#            input_data['Fz'+self.node_id[i]] = print_float_8(f_node[i, 2])
            
        #Assign each node coordiantes to its corresponding node ID in the input data dictionary
        for i in range(len(node_coord_all)):
            input_data['x'+self.node_id_all[i]] = print_float_8(node_coord_all[i,0])
            input_data['y'+self.node_id_all[i]] = print_float_8(node_coord_all[i,1])
            input_data['z'+self.node_id_all[i]] = print_float_8(node_coord_all[i,2])
            
        #Assign each thickness value to its corresponding ID in the input data dictionary
        for i in range(len(t)):
            input_data['t'+str(i+1)] = print_float_8(t[i])
            
        #Assign each stringer section value to its corresponding ID in the input data dictionary
        for i in range(len(s)):
            input_data['s'+str(i+1)] = print_float_8(s[i])
            
        #Assign the Young's modulus to its input data dictionary key
        input_data['E'] = print_float_8(E)
        
        #Assign the Poisson's ratio to its input data dictionary key
        input_data['nu'] = print_float_8(nu)
        
        #Assign the material density to its input data dictionary key
        input_data['rho_s'] = print_float_8(rho_s)
        
#        #In case the stiffnes matrix has not changed, load the stored stiffness matrix. Otherwise, compute and store the new one
#        if not changed:
#            input_data['K'] = 'ALTER 82,82\nINPUTT2   /LLL,,,,/C,N,-1/C,N,11\nENDALTER'
#        else:
#            input_data['K'] = 'ALTER 82\nOUTPUT2  LLL,,,,//C,N,-1 / C,N,11\nENDALTER'
        
                
        #Read the input file template
        f = open(self.template_file,'r')
        tmp = f.read()
        f.close()

        #Replace the input data contained in the dictionary onto the new input file       
        new_file = tmp.format(**input_data)

        inp = open(self.input_filepath,'w')
        inp.write(new_file)
        inp.close()
        
        
#    def get_output_data(self):
#        
#        #Read the punch and output files only if they exist and their last modified date is older than input file one
#        
#        while(not os.path.isfile(self.output_filepath)): pass
#        
#        while(os.path.getmtime(self.output_filepath) <= os.path.getmtime(self.input_filepath)): pass
#    
#        while(not os.path.isfile(self.output_file)): pass
#        
#        while(os.path.getmtime(self.output_file) <= os.path.getmtime(self.input_filepath)): pass
#        
#        u = np.zeros((self.ns,3))
#        
#        shell_stress = []
#        
#        rod_stress = []
#        
#        #Read the Nastran punch file (.pnh) and extract displacement and stress data        
#        with open(self.output_filepath) as f:
#            lines = f.readlines()
#            lines = [i.split() for i in lines]
#            
#            rod_type = False
#            for i in range(len(lines)):
#                if len(lines[i]) > 1:
#                    #Set rod_type to True if the begining of the rod type section is reached
#                    if len(lines[i]) > 3:
#                        if lines[i][0] == '$ELEMENT' and lines[i][4] == '(ROD':
#                            #print(lines[i][4])
#                            rod_type = True
##                        else:
##                            rod_type = False
#                            
##                    if len(lines[i]) > 2:
##                        if lines[i][4] == '(ROD':
##                            rod_type = True
##                        elif lines[i][0] == '$ELEMENT TYPE':
##                            rod_type = False          
#                            
#                    #Write nodal displacements onto u if the node belongs to the outer surface
#                    if lines[i][0] in self.node_id and lines[i][1] == 'G':
#                        u[self.node_id.index(lines[i][0])][0] = lines[i][2]
#                        u[self.node_id.index(lines[i][0])][1] = lines[i][3]
#                        u[self.node_id.index(lines[i][0])][2] = lines[i][4]
#                        
#                    if isint(lines[i][0]) and isfloat(lines[i][1]):   
#                        if rod_type:
#                            #Write rod Von Mises stress onto a list
#                            rod_stress.append(float(lines[i][1]))
#                    
#                        else:
#                            #Write shell principal stresses onto a list (upper and lower shell faces)
#                            shell_stress.append(((float(lines[i+1][3]), float(lines[i+2][1])), (float(lines[i+4][2]), float(lines[i+4][3]))))
#                            
#        #Compute the Von Mises Stress on the structure                
#        VM = []
#                
#        for s in shell_stress:
#            VM.append(np.sqrt(s[0][0]**2 - s[0][0]*s[0][1] + s[0][1]**2))
#            VM.append(np.sqrt(s[1][0]**2 - s[1][0]*s[1][1] + s[1][1]**2))
#            #VM.append(np.sqrt(s[0]**2 - s[0]*s[1] + s[1]**2))
#            
#        for s in rod_stress:
#            VM.append(abs(s))
#        
#        VMStress = np.asarray(VM)        
#        
#        #Read the Nastran output file (.out) and extract the total mass of the structure (M)
#        with open(self.output_file) as f:
#            lines = f.readlines()
#            lines = [i.split() for i in lines]
#            
#            for i in range(len(lines)):
#                if len(lines[i]) > 4:
#                    if lines[i][4] == 'MASS' and lines[i][5] == 'X-C.G.':
#                        M = float(lines[i+1][1].replace('D', 'E'))
#                        
#        output_data = {}                
#        
#        output_data['u'] = u
#        output_data['VMStress'] = VMStress
#        output_data['M'] = M
#        
#        return output_data