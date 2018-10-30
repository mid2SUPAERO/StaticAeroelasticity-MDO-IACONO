# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 10:39:15 2017

@author: c.bruni, j.mascolomer
"""

from __future__ import print_function

#from openmdao.api import ExternalCode

import numpy as np

import os.path, time

from primitives_struct import StructAirfoil

#from primitives_aero import AeroAirfoil

#from liftingsurface import LiftingSurface

from structuralelements import StructuralElements

from OCC.GeomAbs import GeomAbs_C0

import aocxchange.iges

import re

import os

class Geometry(object):
    struct_igs_file = 'struct.igs'
    
    aero_igs_file = 'aero.igs'
    
    struct_geo_file = 'geom.geo'
    
    struct_mesh_file = 'struct_mesh.bdf'
    
    
    def __init__(self, scale_factor, chord_factor ):
        #super(Geometry, self).__init__()
        
        #Number of ribs
        self.n_ribs = 21
        
        #Number of wing sections
        self.sec_num = self.n_ribs-1
        
        #Total number of structure nodes
        self.ns_all = 987
        
        #ID of all nodes
        self.node_id_all = map(str,range(1,self.ns_all+1))
                
        #Number of leading edge stiffeners
        self.n_stiff_le = 3
        
        #Number of trailing edge stiffeners
        self.n_stiff_te = 2              
        
        #Vector containing twist of each section
        self.theta = np.zeros(self.sec_num)
        
        #Global scale factor
        self.scale_factor = scale_factor
        
        #Chord scale factor
        self.chord_factor = chord_factor
        
        #Wing sweep angle
        self.sweep = 0.
        
        #Mid spar position
        self.mid_spar_pos = 0.5
        
        #Coordinates of all structural nodes
        self.node_coord_all = np.zeros((self.ns_all, 3))
        
        self.input_filepath  = self.struct_geo_file
        self.output_filepath = self.struct_mesh_file
        
        #Check if the files exist (optional)        
        #self.options['external_input_files'] = [self.input_filepath,]
        
        #self.options['command'] = ['cmd.exe', '/c', r'gmsh.exe', self.input_filepath, r'-0', r'-o', self.output_filepath]
        
        
    def solve_nonlinear(self):
        
        # Generate the geometry (.stp and .igs files) from the input parameters
        self.create_geom()
        
        # Generate the input file for GMSH from the input parameters
        self.create_struct_mesh()
        
        # Parent solve_nonlinear function actually runs the external code
        #super(Geometry, self).solve_nonlinear()
        self.run_gmsh()
#        
#        #Read aerodynamic points coordinates
#        aero_output_data = self.get_aero_output_data()
#        
        #Read structure points coordinates
        struct_output_data = self.get_struct_output_data()
#        
#        #Write aerodynamic points into output dictionary
#        unknowns['apoints_coord'] = aero_output_data['apoints_coord']
#        
#        #Write outer surface structure points into output dictionary
#        unknowns['node_coord'] = struct_output_data['node_coord']
#        
#        #Write structure points into output dictionary
#        unknowns['node_coord_all'] = struct_output_data['node_coord_all'] 


        return struct_output_data
        
    def create_geom(self):
        
        # Generate the geometry (.stp) from the input parameters
#        theta = params['theta']
#        scale_factor = params['scale_factor']
#        chord_factor = params['chord_factor']
#        sweep = params['sweep']
#        mid_spar_pos = params['mid_spar_pos']
        theta        = self.theta
        scale_factor = self.scale_factor
        chord_factor = self.chord_factor
        sweep        = self.sweep
        mid_spar_pos = self.mid_spar_pos
        
        def myDihedralFunctionAirliner(Epsilon):
            """User-defined function describing the variation of dihedral as a function
            of the leading edge coordinate"""
            BaseDihedral = 0  
            return BaseDihedral 
    
        def myTwistFunctionAirliner(Epsilon):
            """User-defined function describing the variation of sweep angle as a function
            of the leading edge coordinate"""
            Twist = np.zeros(self.sec_num+1)
            #Enforce the root section to zero twist
            Twist[0] = 0.
            Twist[1:] = theta
            EpsArray = np.linspace(0, 1, self.sec_num+1)
                
            return np.interp(Epsilon, EpsArray, Twist)
        
        
        def myChordFunctionAirliner_struct(Epsilon):
            """User-defined function describing the variation of chord as a function of
            the leading edge coordinate"""
            ChordLengths = np.array([0.324466566, 0.28351649, 0.263040242, 0.242562784, 0.222085325, 0.201607867, 0.181129803, 0.172939183, 0.167799092, 0.159231467, 0.150664448, 0.142097428,
              0.133529803, 0.124962784, 0.116395764, 0.107828139, 0.09926112, 0.0906941,0.082126475, 0.073559455, 0.064992436])
            EpsArray = np.linspace(0, 1, 21)
            return np.interp(Epsilon, EpsArray, ChordLengths)
            
            
        def myChordFunctionAirliner_aero(Epsilon):
            """User-defined function describing the variation of chord as a function of
            the leading edge coordinate"""
        
            ChordLengths = np.array([0.463, 0.405, 0.375771774, 0.346518262, 0.31726475, 0.288011238, 0.258756862, 0.247055976, 0.239712989, 0.227473525, 0.215234925,
              0.202996326, 0.190756862, 0.178518262, 0.166279663, 0.154040199, 0.141801599, 0.129563, 0.117323536, 0.105084936, 0.092846337])
            EpsArray = np.linspace(0, 1, 21)
            return np.interp(Epsilon, EpsArray, ChordLengths)
        
        def myAirfoilFunctionAirliner_struct(Epsilon, LEPoint, ChordFunct, ChordFactor,
                                      DihedralFunct, TwistFunct):
            """Defines the variation of cross section as a function of Epsilon"""
            AfChord = ((ChordFactor*ChordFunct(Epsilon)) /
                       np.cos(np.radians(TwistFunct(Epsilon))))
                   
            Af = StructAirfoil(LEPoint, ChordLength=AfChord,
                                    Rotation=DihedralFunct(Epsilon),
                                    Twist=TwistFunct(Epsilon),
                                    CRMProfile=True, CRM_Epsilon=Epsilon)
                                      
            return Af
        
        
        def myAirfoilFunctionAirliner_aero(Epsilon, LEPoint, ChordFunct, ChordFactor,
                                      DihedralFunct, TwistFunct):
            """Defines the variation of cross section as a function of Epsilon"""
            AfChord = ((ChordFactor*ChordFunct(Epsilon)) /
                       np.cos(np.radians(TwistFunct(Epsilon))))
                   
            Af = AeroAirfoil(LEPoint, ChordLength=AfChord,
                                    Rotation=DihedralFunct(Epsilon),
                                    Twist=TwistFunct(Epsilon),
                                    CRMProfile=True, CRM_Epsilon=Epsilon)
                                      
            return Af    
        
        def mySweepAngleFunctionAirliner_struct(Epsilon):
            """User-defined function describing the variation of sweep angle as a function
            of the leading edge coordinate"""
            SweepAngles = sweep*np.ones(self.sec_num+1)
            EpsArray = np.linspace(0, 1, self.sec_num+1)
                
            return np.interp(Epsilon, EpsArray, SweepAngles)
            
        
        def mySweepAngleFunctionAirliner_aero(Epsilon):
            """User-defined function describing the variation of sweep angle as a function
            of the leading edge coordinate"""
            SweepAngles = sweep*np.ones(self.sec_num+1)
            EpsArray = np.linspace(0, 1, self.sec_num+1)
                
            return np.interp(Epsilon, EpsArray, SweepAngles)
       
        # Position of the apex of the wing
        P = (0,0,0)
        NRibs=self.n_ribs
        Multiple=1
        NSegLoft = Multiple*(NRibs-1)
        PosLESpar=0.0 #Chord percentage
        PosMIDSpar=mid_spar_pos
        PosTESpar=1.0 #Chord percentage
        NStiffLE=self.n_stiff_le
        NStiffTE=self.n_stiff_te
        ChordFactor = chord_factor
        ScaleFactor = scale_factor
#        NP_chord = 2*self.network_info[0][1] - 1
#        NP_span = self.network_info[0][2] - 1
        
        
#        Wing = LiftingSurface(P, mySweepAngleFunctionAirliner_aero, 
#                                     myDihedralFunctionAirliner, 
#                                     myTwistFunctionAirliner, 
#                                     myChordFunctionAirliner_aero, 
#                                     myAirfoilFunctionAirliner_aero, 
#                                     SegmentNo=NSegLoft,
#                                     NPaero_chord=NP_chord,
#                                     NPaero_span=NP_span,
#                                     ChordFactor=ChordFactor,
#                                     ScaleFactor=ScaleFactor)
                                     
                                  
        WB=StructuralElements(P, mySweepAngleFunctionAirliner_struct, 
                                     myDihedralFunctionAirliner, 
                                     myTwistFunctionAirliner, 
                                     myChordFunctionAirliner_struct, 
                                     myAirfoilFunctionAirliner_struct, 
                                     SegmentNoLoft=NSegLoft,
                                     NoRibs=NRibs,
                                     PositionLESpar=PosLESpar,
                                     PositionMIDSpar=PosMIDSpar,
                                     PositionTESpar=PosTESpar,
                                     NoStiffnersLE=NStiffLE,
                                     NoStiffnersTE=NStiffTE,
                                     ChordFactor=ChordFactor,
                                     ScaleFactor=ScaleFactor,
                                     max_degree=8,
                                     continuity=GeomAbs_C0,
                                     )
                                     
        #Get current working directory
        cwd = os.getcwd()
        
#        # Export IGES Aero
#        filename = aocxchange.utils.path_from_file(__file__, cwd + '//' + self.aero_igs_file)
#        my_iges_exporter = aocxchange.iges.IgesExporter(filename, format="5.1")
#        PointsUp=Wing['PointsSurfUp']
#        PointsDown=Wing['PointsSurfDown']
#        my_iges_exporter.add_shape(PointsUp)
#        my_iges_exporter.add_shape(PointsDown)   
#        my_iges_exporter.write_file() 
        
        # Export IGES Struct
        filename = aocxchange.utils.path_from_file(__file__, cwd + '//' + self.struct_igs_file)
        my_iges_exporter = aocxchange.iges.IgesExporter(filename, format="5.1")
        Ribs=WB['Ribs']   
        my_iges_exporter.add_shape(Ribs)
        Panels=WB['Panels']
        my_iges_exporter.add_shape(Panels)    
        Spars=WB['Spars']
        my_iges_exporter.add_shape(Spars)
        Stiffners=WB['Stringers']
        my_iges_exporter.add_shape(Stiffners)
        my_iges_exporter.write_file()
            
      
    def create_struct_mesh(self):
        
        #Open struct .igs file and get total number of surfaces
        with open(self.struct_igs_file) as f:
            lines = f.readlines()
            lines = [i.split() for i in lines]
            
        #Surface counter
        N_surf = 0
        
        #Count double of number of surfaces
        for line in lines:
            if line[0] == '144':
                N_surf += 1;
                
        #Write .geo file
        with open(self.struct_geo_file,'w') as f:
            f.write('Merge "'+self.struct_igs_file+'";\n')
            f.write('Transfinite Surface {')
            for i in range(1, N_surf/2):
                f.write(str(i)+',')
            f.write(str(N_surf/2)+'};\n')
            f.write('Recombine Surface {')
            for i in range(1, N_surf/2):
                f.write(str(i)+',')
            f.write(str(N_surf/2)+'};\n')
            f.write('Mesh.RecombineAll = 1;')
            f.write('Mesh 2;  // Generate 2D mesh\n')
            f.write('Coherence Mesh;  // Remove duplicate entities\n')
            f.write('Mesh.Format = 31;\n')
            f.write('Mesh.BdfFieldFormat = 0;\n')             
            f.write('Save "'+self.struct_mesh_file+'";\n')
   


    def run_gmsh(self):
        os.system('gmsh '+self.struct_geo_file+' -0')


         
            
    def get_struct_output_data(self):
        
        struct_output_data = {}
        
#        node_coord = np.zeros((len(self.node_id),3))
        node_coord_all = np.zeros((len(self.node_id_all),3))
        
        #Read total and outer surface node coordinates from file and write them into an array        
        with open(self.struct_mesh_file) as f:
            lines = f.readlines()
            lines = [re.split(r'[,\s]\s*', i) for i in lines]
        
        for line in lines:
            if len(line) > 1:
                if line[0] == 'GRID' and line[1] in self.node_id_all:
                    node_coord_all[self.node_id_all.index(line[1]), 0] = float(line[3])                    
                    node_coord_all[self.node_id_all.index(line[1]), 1] = float(line[4])
                    node_coord_all[self.node_id_all.index(line[1]), 2] = float(line[5])
#                    if line[1] in self.node_id:
#                        node_coord[self.node_id.index(line[1]), 0] = float(line[3])
#                        node_coord[self.node_id.index(line[1]), 1] = float(line[4])
#                        node_coord[self.node_id.index(line[1]), 2] = float(line[5])
        
#        struct_output_data['node_coord'] = node_coord
        struct_output_data['node_coord_all'] = node_coord_all
        
        return struct_output_data
        
        
#    def get_aero_output_data(self):
#        
#        aero_output_data = {}
#        
#        #Read aero points coordinates
#        with open(self.aero_igs_file) as f:
#            lines = f.readlines()
#            lines = [i.split(',') for i in lines]
#        
#        apoints_coord = []
#        
#        #Read aero points coordinates.
#        for line in lines:
#            if len(line) > 4:
#                if line[0] == '116' and line[4][:2] == '0;':
#                    apoints_coord.append([float(line[1]), float(line[2]), float(line[3])])
#
#        apoints_coord = np.asarray(apoints_coord)
#        
#        aero_output_data['apoints_coord'] = apoints_coord
#        
#        return aero_output_data
