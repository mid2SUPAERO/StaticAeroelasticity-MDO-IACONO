# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 08:59:43 2016

@author: jmascolo
"""

from __future__ import print_function

import numpy as np

#from is_number import isfloat, isint

class ProblemData:
    #RBF function type for interpolation
    function_type = 'thin_plate'
    
    #Epsilon parameter of the RBF functions
    epsilon = None
    
    #Norm bias
    bias = None
    
    #Jig shape geometry
    aero_template = 'aero_template.wgs'
    
    #Nastran template file
    template_file = 'nastran_input_template.inp'
    
    #Number of wing sections
    sec_num = 20
    
    #Number of ribs
    n_ribs = 21
    
    #Number of leading edge stiffeners
    n_stiff_le = 3
    
    #Number of trailing edge stiffeners
    n_stiff_te = 2
    
    
    def __init__(self):
        #Number of aerodynamic grid points
        self.na = self.get_aero()['na']
        
        #Aerodynamic mesh network information list
        self.network_info = self.get_aero()['network_info']
        
        #List of structural nodes identification numbers (on the outer surface)
        self.node_id = self.get_structure()['node_id']
        
        #List of total structural nodes identification numbers (on the outer surface)
        self.node_id_all = self.get_structure()['node_id_all']
        
        #Number of different thickness regions
        self.tn = self.get_structure()['tn']
        
        #Number of different stringer sections
        self.sn = self.get_structure()['sn']
        
        #Number of structural nodes on the outer surface
        self.ns = len(self.node_id)
        
        #Total number of structural nodes
        self.ns_all = len(self.node_id_all)
        
        #Number of Von Mises stress outputs
        self.n_stress = self.get_structure()['n_stress']
        
        
    def get_aero(self):
        
        apoints_coord = []
        network_info = []
        
        #Write the aerodynamic grid points coordinates into a list
        with open(self.aero_template) as f:
            lines = f.readlines()
            lines = [i.split() for i in lines]
            
            points = 0
            pans = 0
            
            for line in lines:
                #Get aerodynamic grid points coordinates
                if all(isfloat(item) for item in line):
                    if len(line) == 3:
                        apoints_coord.append([float(line[0]), float(line[1]), float(line[2])])
                    if len(line) == 6:
                        apoints_coord.append([float(line[0]), float(line[1]), float(line[2])])
                        apoints_coord.append([float(line[3]), float(line[4]), float(line[5])])
                            
                #Get network information
                if len(line) > 1:
                    if isint(line[0]):
                        # I have put it at the end because I think otherwise it cannot get the good info from upper surface
                        network_info.append([int(line[0]), int(line[2]), int(line[1]), int(points), int(pans)])
                        points = points + int(line[2])*int(line[1])
                        pans = pans + (int(line[2]) - 1)*(int(line[1]) - 1)
                        #network_info.append([int(line[0]), int(line[2]), int(line[1]), int(points), int(pans)])
        apoints_coord = np.asarray(apoints_coord)
        na = len(apoints_coord)
        
        #Dictionary containing aerodynamic problem data
        aero = {}
        aero['na'] = na
        aero['network_info'] = network_info
        
        return aero
        
        
    #Function that returns the list of node IDs belonging to the outer surface    
    def get_structure(self):
        
        node_id = []
        node_id_all = []
        tn = 0
        sn = 0
        n_stress = 0
        
        #Read the list of nodes belonging to the outer surface from the template file
        with open(self.template_file) as f:
            lines = f.readlines()
            
            outer_node_begin = lines.index('$List of nodes belonging to the outer skin\n')
            total_node_begin = lines.index('$List of total nodes\n')
            
            for i in range(len(lines)):
                #Store nodes belonging to the outer skin
                if i > outer_node_begin and lines[i][0] == '$':
                    node_id.append(lines[i].strip().lstrip('$'))
                #Store all structure nodes
                if i > total_node_begin and i < outer_node_begin and lines[i][0] == '$':
                    node_id_all.append(lines[i].strip().lstrip('$'))
                #Store number of different thickness regions
                lines[i] = lines[i].split(',')
                if lines[i][0] == 'PSHELL':
                    tn += 1
                #Store number of different stringer sections
                elif lines[i][0] == 'PROD':
                    sn += 1
                #Store number of stress outputs
                elif lines[i][0] == 'CTRIA3' or lines[i][0] == 'CQUAD4':
                    n_stress += 2
                elif lines[i][0] == 'CROD':
                    n_stress += 1
                    
        #Dictionary containing structure problem data
        structure = {}
        structure['node_id'] = node_id
        structure['node_id_all'] = node_id_all
        structure['tn'] = tn
        structure['sn'] = sn
        structure['n_stress'] = n_stress
                    
        return structure