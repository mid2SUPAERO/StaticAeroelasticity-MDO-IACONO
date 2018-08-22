# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 17:20:42 2016 

@author: © Joan Mas Colomer
"""

from __future__ import print_function

import numpy as np

from openmdao.api import Component

'''
Component which takes the nodal displacements and gives the displacements of
the aerodynamic points
'''
class DisplacementTransfer(Component):


    def __init__(self, na, ns):
        super(DisplacementTransfer, self).__init__()

        #Number of points of the aerodynamic grid
        self.na = na

        #Number of nodes of the structural mesh on the outer skin
        self.ns = ns

        #Interpolation matrix H (xa = H xs)
        self.add_param('H', val=np.zeros((self.na, self.ns)))

        #Nodal displacements of the outer surface
        self.add_param('u', val=np.zeros((self.ns, 3)))

        #Displacements of the aerodynamic grid points
        self.add_output('delta', val=np.zeros((self.na, 3)))


    def solve_nonlinear(self, params, unknowns, resids):

        u = params['u']

        H = params['H']

        #Apply the interpolation matrix to obtain the aerodynamic points displacements
        unknowns['delta'] = H.dot(u)
        
#==============================================================================
# class Aggregation(Component):
#     
#     def __init__(self,n_stress,p,StressAllow):
#         super(Aggregation,self).__init__()
#         self.n_stress=n_stress
#         self.p=p
#         self.StressAllow=StressAllow
#         self.add_param('VMStress',np.zeros(self.n_stress))
#         self.add_param('gi',np.zeros(self.n_stress))
#         self.add_param('gmax',val=0.0)
#         self.add_param('summ',val=0.0)
#         self.add_output('Gksl',val=0.0)
#     def solve_nonlinear(self,params,unknowns,resids):
#         
#         n_stress=self.n_stress
#         p=self.p
#         gi=params['gi']
#         VMStress=params['VMStress']
#         gmax=params['gmax']
#         gi=(VMStress/self.StressAllow -1)
#         gmax=max(gi)
#         #print('gmax = '+str(gmax))
#         summ=params['summ']
#         for k in range (n_stress):
#             
#             summ=summ+np.exp(p*(gi[k]-gmax))
#             
#                 
#         unknowns['Gksl']=gmax+(np.log(summ)/p)-(np.log(n_stress)/p)
#         #print ('Gksl =  '+str(unknowns['Gksl']))
#         
# 
#==============================================================================
