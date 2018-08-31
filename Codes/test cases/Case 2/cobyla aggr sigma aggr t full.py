# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 10:50:10 2016

@author: © a.iacono
"""

from __future__ import print_function

from openmdao.api import Problem, Group, IndepVarComp, ScipyGMRES, SqliteRecorder, ExecComp, ScipyOptimizer, view_tree

from aerostructures import NastranStatic, DisplacementTransfer, Panair, LoadTransfer,Aggregation, Interpolation, StaticStructureProblemDimensions, StaticStructureProblemParams, AeroProblemDimensions, AeroProblemParams, NLGaussSeidel

from ttt import Ttt
import numpy as np

if __name__ == "__main__":

    #Interpolation function type and setup
    function_type = 'thin_plate'
    bias = (1,50,1)

    #Symmetry plane index
    sym_plane_index = 1

    #Problem parameters
    Sw = 383.689555
    V = 250.75
    rho_a = 0.337
    Mach = 0.85
#    alpha = 0.58465
    b = 58.7629
    c = 7.00532
    E = 6.89e10
    nu = 0.31
    rho_s = 2795.67
    sigma_y = 5.033172e+08
    t_ip=0.01
    W = 9.81*300000.
    function = 'Gksl'
    p=100.
    s0=40000000.0
    s1=0.01

    structure_problem_dimensions = StaticStructureProblemDimensions()
    aero_problem_dimensions = AeroProblemDimensions()

    ns = structure_problem_dimensions.ns
    ns_all = structure_problem_dimensions.ns_all
    node_id = structure_problem_dimensions.node_id
    node_id_all = structure_problem_dimensions.node_id_all
    n_stress = structure_problem_dimensions.n_stress
    n_t=12
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

    top = Problem()
    top.root = root = Group()
#    top.root.deriv_options['type'] = 'fd'
#    top.root.deriv_options['step_size'] = 1.0e-1
    #Add independent variables
    root.add('wing_area', IndepVarComp('Sw', Sw), promotes=['*'])
    root.add('airspeed', IndepVarComp('V', V), promotes=['*'])
    root.add('sigma_y', IndepVarComp('sigma_y', sigma_y), promotes=['*'])
    root.add('t_ip', IndepVarComp('t_ip', t_ip), promotes=['*'])
    root.add('stress_ref', IndepVarComp('s0', s0), promotes=['*'])
    root.add('tick_ref', IndepVarComp('s1', s1), promotes=['*'])
    root.add('air_density', IndepVarComp('rho_a', rho_a), promotes=['*'])
    root.add('Mach_number', IndepVarComp('Mach', Mach), promotes=['*'])
    root.add('young_module', IndepVarComp('E', E), promotes=['*'])
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

    root.add('inter', Interpolation(na, ns, function = function_type, bias = bias), promotes=['*'])
    root.add('agrr', Aggregation(n_stress,p,function), promotes=['*'])
    root.add('agrrt', Ttt(n_t,p,function), promotes=['*'])
    
    mda = Group()

    #Add disciplines to the group
    mda.add('displacement_transfer', DisplacementTransfer(na, ns), promotes=['*'])
    mda.add('aerodynamics', Panair(na, network_info), promotes=['*'])
    mda.add('load_transfer', LoadTransfer(na, ns), promotes=['*'])
    mda.add('structures', NastranStatic(node_id, node_id_all, n_stress, tn, mn), promotes=['*'])

    #Define solver type and tolerance for MDA
    mda.nl_solver = NLGaussSeidel()
#    mda.nl_solver.options['rtol'] = 1.e-1
    mda.nl_solver.options['maxiter'] = 15
    mda.nl_solver.options['rutol'] = 1.e-2
    mda.nl_solver.options['use_aitken'] = True
    mda.nl_solver.options['aitken_alpha_min'] = 0.1
    mda.nl_solver.options['aitken_alpha_max'] = 1.5

    mda.ln_solver = ScipyGMRES()

    root.add('mda_group', mda, promotes=['*'])

#    top.root.mda_group.deriv_options['type'] = 'fd'
#    top.root.mda_group.deriv_options['step_size'] = 1.0e-1
    #Recorder
    recorder = SqliteRecorder('opti_g_46')
    recorder.options['record_params'] = False
    recorder.options['record_metadata'] = False
    recorder.options['record_resids'] = False
    recorder.options['record_derivs'] = False
    top.root.nl_solver.add_recorder(recorder)
#    top.root.add_recorder(recorder)


    #Define solver type

    root.add('obj_function', ExecComp('obj_f = mass'), promotes=['*'])
    
    t_max=0.01*np.ones(tn)
    t_min=0.006*np.ones(tn)

    
    #Add constraint components
    root.add('con_lift', ExecComp('con_l = CL - W/(0.5*rho_a*V**2*Sw)'), promotes=['*'])
    root.add('con_stress', ExecComp('con_s = G - sigma_y'), promotes=['*'])

    root.add('con_t', ExecComp('con_t = t_ip - Gt'),  promotes=['*'])
    #Define solver type
    root.ln_solver = ScipyGMRES()
    
    #Define the optimizer
    top.driver = ScipyOptimizer()
    top.driver.options['optimizer'] = 'COBYLA'
    top.driver.options['disp'] = True
    top.driver.options['tol'] = 1.e-4
#    top.driver.options['maxiter'] = 500
    
    top.driver.opt_settings['rhobeg']= 0.025
#    top.driver.eps = 0.1
    
    alpha_max=10.
    alpha_min=0.
    

    top.driver.add_desvar('alpha', lower=alpha_min, upper=alpha_max, adder=-alpha_min, scaler=1/(alpha_max-alpha_min))
    top.driver.add_desvar('t', lower=t_min, upper=t_max, adder=-t_min, scaler=1/(t_max-t_min))
    
#    top.driver.add_desvar('alpha', lower=alpha_min, upper=alpha_max)
#    top.driver.add_desvar('t', lower=t_min, upper=t_max)
    top.driver.add_objective('obj_f')
    top.driver.add_constraint('con_l', lower=0.) #CL>0.72
#    top.driver.add_constraint('con_s', upper=0.) #sigma_max<sigma_y

    top.driver.add_constraint('con_t', lower=0.)
    
    top.setup()
    
    view_tree(top, show_browser=False)
    
    top['alpha'] = 2.0
#    top['t'] = 0.0075*np.ones(tn)
    
    top.run()
    top.cleanup()