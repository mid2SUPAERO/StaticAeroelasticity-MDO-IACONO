# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 17:20:42 2016 

@author: © a.iacono
"""

from __future__ import print_function

import numpy as np

from openmdao.api import Component

#Component which gives the aggregation function (G) given the Von Mises Stresses
class Aggregation(Component):
    
    
    def __init__(self,n_stress,p,function):
        super(Aggregation,self).__init__()
        
        #Number of Von Mises Stresses obtained
        self.n_stress=n_stress
        
        #Draw-down factor
        self.p=p
        
        #Aggregation function type
        self.function=function
        
        #Von Mises stresses vector
        self.add_param('VMStress',np.zeros(self.n_stress))
                
        #Reference vqlue for normalization of stresses
        self.add_param('s0',val=0.0)
        
        #Aggregation function
        self.add_output('G',val=0.0)
        
        
    def solve_nonlinear(self,params,unknowns,resids):
        
        n_stress=self.n_stress
        
        p=self.p
        
        function=self.function
        
        VMStress=params['VMStress']
        
        s0=params['s0']
        
        gmax=max(VMStress)
        
        summ=0.0
        
        
        if function=='Gksl':
        
            for k in range (n_stress):
                
                summ+=np.exp(p*((VMStress[k]-gmax)/s0))
                
            G=((gmax/s0)+(np.log(summ)/p)-(np.log(n_stress)/p))*s0
                    
        
        elif function=='Gksu':
            
            for k in range (n_stress):
                
                summ+=np.exp(p*((VMStress[k]-gmax)/s0))
                
                G=((gmax/s0)+(np.log(summ)/p))*s0
                
                
        elif function=='Gpn':
            
            for k in range(n_stress):
                
                summ+=(VMStress[k]/gmax)**p
                
            G=gmax*(summ)**(1/p)   
            
                
        elif function=='Gpm':
            
            for k in range(n_stress):
                
                summ+=(VMStress[k]/gmax)**p
                
            G=(summ/n_stress)**(1/p)*gmax  

               
        #Set the aggregation function as an output
        unknowns['G']=G
        
        output_aggr = {}

        output_aggr['G'] = G


        return output_aggr


                

