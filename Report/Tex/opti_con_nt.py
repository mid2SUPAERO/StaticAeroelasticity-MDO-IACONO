#MAIN SCRIPT
#Modules Importing
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 10:50:10 2016

@author:  a.iacono
"""
from __future__ import print_function

from openmdao.api import Problem, Group, IndepVarComp, ScipyGMRES, SqliteRecorder, ExecComp, ScipyOptimizer, view_tree

from aerostructures import NastranStatic, DisplacementTransfer, Panair, LoadTransfer,Aggregation, Interpolation,\
StaticStructureProblemDimensions, StaticStructureProblemParams, AeroProblemDimensions, AeroProblemParams, NLGaussSeidel

import numpy as np
#Definition of Specific Problem Parameters
#Note: To choose correctly the aggregation function parameters checking the guide.

if __name__ == "__main__":

    #Interpolation function type and setup
    function_type = 'thin_plate'
    bias = (1,50,1)

    #Symmetry plane index
    sym_plane_index = 1

    #Problem parameters
    Sw = 383.689555        # Wing Surface
    V = 250.75             # Velocity
    rho_a = 0.337          # Air Density
    Mach = 0.85            # Flight MACH
    alpha = 0.58465        # Starting Angle of Attack
    b = 58.7629            # Wing Span
    c = 7.00532            # Cord
    E = 6.89e10            # Young Module
    nu = 0.31              # Poisson Coefficient
    rho_s = 2795.67        # Material Density
    t_i_max=0.0125         # Upper Limit for Panel's Thickness
    t_i_min=0.0018         # Lower Limit for Panel's Thickness
    sigma_y = 5.033172e+08 # Yield Stress
    W = 9.81*300000.       # Airplane Weight
    function = 'Gksl'      # Aggregation Function
    p=100.                 # Draw-Down Function
    s0=40000000.0          # Reference Stress
#Definition of Problem Dimensions:
#Note: To do this we are using some of the methods defined in the Aerostructures Package, in the relative guide is declared the meaning of each variables and how we determinate it.

    structure_problem_dimensions = StaticStructureProblemDimensions()
    aero_problem_dimensions = AeroProblemDimensions()

    ns = structure_problem_dimensions.ns
    ns_all = structure_problem_dimensions.ns_all
    node_id = structure_problem_dimensions.node_id
    node_id_all = structure_problem_dimensions.node_id_all
    n_stress = structure_problem_dimensions.n_stress
    tn = structure_problem_dimensions.tn
    mn = structure_problem_dimensions.mn

    structure_problem_params = StaticStructureProblemParams(node_id, node_id_all)
    aero_problem_params = AeroProblemParams()

    na = aero_problem_dimensions.na
    network_info = aero_problem_dimensions.network_info

    node_coord = structure_problem_params.node_coord
    node_coord_all = structure_problem_params.node_coord_all
    t = structure_problem_params.t
    m = structure_problem_params.m

    apoints_coord = aero_problem_params.apoints_coord
#Setting OpenMDAO Problem

    top = Problem()
    top.root = root = Group()
#If you want to used a gradient based optimization like SLSQP the gradient is obtained by finite difference. So you have to set the step because the default step make really little changes on the variables, and the 8 float architecture of nastran can't recognize this variation.

#=================================================
#     UNCOMMENT JUST FOR SLSQP OPTIMIZER
#=================================================
#     top.root.deriv_options['type'] = 'fd'
#     top.root.deriv_options['step_size'] = 1.0e-1
#=================================================
#Here we are defining the openMDAO variables as Independent Variables Component:

    root.add('wing_area', IndepVarComp('Sw', Sw), promotes=['*'])
    root.add('airspeed', IndepVarComp('V', V), promotes=['*'])
    root.add('sigma_y', IndepVarComp('sigma_y', sigma_y), promotes=['*'])
    root.add('stress_ref', IndepVarComp('s0', s0), promotes=['*'])
    root.add('air_density', IndepVarComp('rho_a', rho_a), promotes=['*'])
    root.add('Mach_number', IndepVarComp('Mach', Mach), promotes=['*'])
    root.add('young_module', IndepVarComp('E', E), promotes=['*'])
    root.add('tick_max',IndepVarComp('t_i_max',t_i_max),promotes=['*'])
    root.add('t_min', IndepVarComp('t_i_min',t_i_min), promotes=['*'])
    root.add('weight', IndepVarComp('W', W), promotes=['*'])
    root.add('mat_density', IndepVarComp('rho_s', rho_s), promotes=['*'])
    root.add('poisson', IndepVarComp('nu', nu), promotes=['*'])
    root.add('angle_of_attack', IndepVarComp('alpha', 0.), promotes=['*'])
    root.add('wing_span', IndepVarComp('b', b), promotes=['*'])
    root.add('wing_chord', IndepVarComp('c', c), promotes=['*'])
    root.add('s_coord', IndepVarComp('node_coord', node_coord), promotes=['*'])
    root.add('s_coord_all', IndepVarComp('node_coord_all', node_coord_all), promotes=['*'])
    root.add('thicknesses', IndepVarComp('t', t), promotes=['*'])
    root.add('masses', IndepVarComp('m', m), promotes=['*'])
    root.add('a_coord', IndepVarComp('apoints_coord', apoints_coord), promotes=['*'])
#Here we are adding the modules for the interpolation between aerodinamic and strucutural mesh, and the aggregation module, defined in the Aereostructures package:

    root.add('inter', Interpolation(na, ns, function = function_type, bias = bias), promotes=['*'])
    root.add('agrr', Aggregation(n_stress,p,function), promotes=['*'])
    
#Creating the mda group for the convergence of aerodinamyc forces and structural displacement:
    mda = Group()

    #Add disciplines to the group
    mda.add('displacement_transfer', DisplacementTransfer(na, ns), promotes=['*'])
    mda.add('aerodynamics', Panair(na, network_info), promotes=['*'])
    mda.add('load_transfer', LoadTransfer(na, ns), promotes=['*'])
    mda.add('structures', NastranStatic(node_id, node_id_all, n_stress, tn, mn), promotes=['*'])
#Setting the mda solver type and settings:
    mda.nl_solver = NLGaussSeidel()
#    mda.nl_solver.options['rtol'] = 1.e-1
    mda.nl_solver.options['maxiter'] = 15
    mda.nl_solver.options['rutol'] = 1.e-2
    mda.nl_solver.options['use_aitken'] = True
    mda.nl_solver.options['aitken_alpha_min'] = 0.1
    mda.nl_solver.options['aitken_alpha_max'] = 1.5
#Adding the mda cycle to the optimization group:

    mda.ln_solver = ScipyGMRES()

    root.add('mda_group', mda, promotes=['*'])
#====================================================
#      UNCOMMENT JUST FOR SLSQP OPTIMIZER
#====================================================
#====================================================
#     top.root.mda_group.deriv_options['type'] = 'fd'
#     top.root.mda_group.deriv_options['step_size'] = 1.0e-1
#====================================================
#Setting Recorder
    recorder = SqliteRecorder('opti_g_55')       # In brackets the name of the database
    recorder.options['record_params'] = False
    recorder.options['record_metadata'] = False
    recorder.options['record_resids'] = False
    recorder.options['record_derivs'] = False
    top.root.nl_solver.add_recorder(recorder)
#Setting the Optimitazion
#Defining solver type and the optimizer
    root.ln_solver = ScipyGMRES()               # Solver type
    top.driver = ScipyOptimizer()               # Adding Optmizer 
    top.driver.options['optimizer'] = 'COBYLA'  # Defining the optimizer and is settings
    top.driver.options['disp'] = True           
    top.driver.options['tol'] = 1.e-4
    top.driver.options['maxiter'] = 500
    top.driver.opt_settings['rhobeg']= 0.1      # USE JUST FOR COBYLA
#Defining objective and constraint functions
    # Minimize the CDi or the mass
    root.add('obj_function', ExecComp('obj_f = CDi'), promotes=['*'])
    # Define constraint CL - W/q
    root.add('con_lift', ExecComp('con_l = CL - W/(0.5*rho_a*V**2*Sw)'), promotes=['*'])
    # Define constraint G - sigma_y
    root.add('con_stress', ExecComp('con_s = G - sigma_y'), promotes=['*'])                    
    t_max=0.01*np.ones(tn)
    t_min=0.006*np.ones(tn)
    # Define one constraint t_max and t_min for each differentthickness section in the model
    for i in range(tn):
        root.add('max_t_'+str(i+1),ExecComp('max_t_'+str(i+1)+' = t['+str(i)+'] - t_i_max',     
                                            t=np.zeros(tn,dtype=float)),promotes=['*'])         
        root.add('min_t_'+str(i+1),ExecComp('min_t_'+str(i+1)+' = t['+str(i)+'] - t_i_min',    
                                            t=np.zeros(tn,dtype=float)),promotes=['*'])
#Setting objective, design variables and constraint
    top.driver.add_objective('obj_f')
    
    alpha_max=10.
    alpha_min=0.
    
    top.driver.add_desvar('alpha', lower=alpha_min, upper=alpha_max, adder=-alpha_min, scaler=1/(alpha_max-alpha_min))
    top.driver.add_desvar('t', lower=t_min, upper=t_max, adder=-t_min, scaler=1/(t_max-t_min))
    
    top.driver.add_constraint('con_l', lower=0.) 
    top.driver.add_constraint('con_s', upper=0.) 
    for i in range(tn):
        top.driver.add_constraint('max_t_'+str(i+1),upper=0.,scaler=1/t_i_max)
        top.driver.add_constraint('min_t_'+str(i+1),lower=0.,scaler=1/t_i_min)
#As u can see we have created one constraint for each different thickness section, but for
#the stress constraint we just created one constraint using the aggregation function, that is
#representative of the max value of the stress for each iteration.

#Setup Problem and Run It
    top.setup()
    
    view_tree(top, show_browser=False)
    
    top.run()
    top.cleanup()
