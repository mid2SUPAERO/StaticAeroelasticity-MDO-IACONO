# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 15:35:16 2016

@author: c.bruni
"""

import numpy as np
from airconics.base import AirconicsShape
import airconics.AirCONICStools as act
from OCC.gp import gp_Pnt, gp_Vec, gp_Pln, gp_Ax1, gp_Dir, gp_OY
from OCC.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_Line, GeomAbs_C1  
from OCC.Geom import Geom_Plane
from OCC.GC import GC_MakeSegment
from OCC.GeomAPI import GeomAPI_IntCS, GeomAPI_PointsToBSplineSurface
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeFace       
from OCC.BRepAlgoAPI import BRepAlgoAPI_Section  
from OCC.GeomPlate import GeomPlate_BuildPlateSurface,GeomPlate_MakeApprox
from OCC.GeomFill import GeomFill_BSplineCurves  
from OCC.BRepFill import BRepFill_CurveConstraint
from OCC.BRepAdaptor import BRepAdaptor_HCurve 
from OCC.BRepTools import BRepTools_WireExplorer
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.TColgp import TColgp_Array2OfPnt
from OCC.GeomFill import GeomFill_StretchStyle, GeomFill_CoonsStyle, GeomFill_CurvedStyle
#from franges import drange, frange

# The current script allow to build the structural elements such as: Ribs, Spars, Panels and Stringers.
# The class "StructuralElements" is organized in such a way: 
# The method "_GenerateLeadingEdge" provide the position of the LE along the dimensionless wingspan.
# The method "_BuildRibs" allows to create the Ribs along the dimensionless wingspan.
# The method "_BuildSpars" create the Spars and the Stringers along the dimensionless wingspan.The are togheter since they are based on the same 
# leading edge coordinates definition.  
# The method "_BuildPanels" generate the panels with respect to the number of stringer and of ribs.             
                
                           
                                
class StructuralElements(AirconicsShape):
    """Class for defining structural shapes

    Parameters
    ----------
    ApexPoint - array, length 3
        Foremost point of the wing (x direction)

    SweepFunct - function
        function defining the leading edge sweep vs epsilon spanwise 
        variable coordinate between 0 and 1 (curvilinear attached
        coordinates)
        
    DihedralFunct - function
        function defining the leading edge dihedral vs epsilon spanwise 
        variable coordinate between 0 and 1 (curvilinear attached)  

    TwistFunc - function
        function defining the sectional twist vs epsilon spanwise 
        variable coordinate between 0 and 1 (curvilinear attached)

    ChordFunct - function
        function defining the leading edge chord vs epsilon spanwise 
        variable coordinate between 0 and 1 (curvilinear attached)    

    ChordFactor - int (default = 1)
        Scaling factor applied in chordwise direction

    ScaleFactor - int (default = 1)
        Scaling factor applied in all directions (uniform)

       
    SegmentNoLoft - int (default = 10)
        Number of segments to sample the spars
        
        
    NoStiffners - int (default = 10)
        Number of stiffners in for each wing section 
        
        
    Attributes
    ----------
    self['Surface'] : TopoDS_Shape
        The generated lifting surface

    See also
    --------
    airconics.primitives.Airfoil
    
    """
    def __init__(self, ApexPoint,
                 SweepFunct,
                 DihedralFunct,
                 TwistFunct,
                 ChordFunct,
                 AirfoilFunct,
                 ChordFactor=1,
                 ScaleFactor=1,
                 OptimizeChordScale=0,
                 LooseSurf=1,
                 SegmentNoLoft=10,
                 NoRibs=10,
                 PositionLESpar=1,
                 PositionMIDSpar=1,
                 PositionTESpar=1,
                 NoStiffnersLE=1,
                 NoStiffnersTE=1,
                 TipRequired=False,
                 max_degree=8,
                 continuity=GeomAbs_C0,):
#        Initialise the components using base class:
        super(StructuralElements, self).__init__(components={},
                                         ApexPoint=gp_Pnt(*ApexPoint),
                                         SweepFunct=SweepFunct,
                                         DihedralFunct=DihedralFunct,
                                         TwistFunct=TwistFunct,
                                         ChordFunct=ChordFunct,
                                         AirfoilFunct=AirfoilFunct,
                                         ChordFactor=ChordFactor,
                                         ScaleFactor=ScaleFactor,
                                         OptimizeChordScale=OptimizeChordScale,
                                         LooseSurf=LooseSurf,
                                         SegmentNoLoft=SegmentNoLoft,
                                         NoRibs=NoRibs,
                                         PositionLESpar=PositionLESpar,
                                         PositionMIDSpar=PositionMIDSpar,
                                         PositionTESpar=PositionTESpar,
                                         NoStiffnersLE=NoStiffnersLE,
                                         NoStiffnersTE=NoStiffnersTE,
                                         TipRequired=TipRequired,
                                         _max_degree=max_degree,
                                         _Cont=continuity)


        self.GenerateStruct(self.ChordFactor, self.ScaleFactor)
        
    def _GenerateLeadingEdge(self):
        """Epsilon coordinate attached to leading edge defines sweep
         Returns airfoil leading edge points
         The generation of a LE spars different from Ribs comes from teh fact 
         that the spars generation comes with the N° of Loft sections, wich is
         generally higher than the number of Ribs sections. If we link the spars
         generation to the Ribs we may have a spar shape which do not follow exactly the loft.
         """

        SegmentLengthLoft = 1. / self.SegmentNoLoft

#       Array of epsilon at segment midpoints (will evaluate curve here)
        Epsilon_midpointsS = np.linspace(SegmentLengthLoft/2., 1-SegmentLengthLoft/2.,
                                        self.SegmentNoLoft)

        Tilt_arrayS = self.DihedralFunct(Epsilon_midpointsS)
        Sweep_arrayS = self.SweepFunct(Epsilon_midpointsS)

        DeltaXsS =  SegmentLengthLoft * np.sin(Sweep_arrayS*(np.pi/180.))
        DeltaYsS = (SegmentLengthLoft * np.cos(Tilt_arrayS*np.pi/180.) *
                                   np.cos(Sweep_arrayS*np.pi/180.))
        DeltaZsS = DeltaYsS*np.tan(Tilt_arrayS*np.pi/180.)

        LEPointsS = np.zeros((self.SegmentNoLoft + 1, 3))

        DeltasS = np.vstack([DeltaXsS, DeltaYsS, DeltaZsS]).T
        LEPointsS[1:, :] = np.cumsum(DeltasS, axis=0)

        return LEPointsS
        

    
    def _BuildSpars(self, ChordFactor, ScaleFactor):
        
        LE = self._GenerateLeadingEdge()
        
        self.NoStiffners=self.NoStiffnersLE+self.NoStiffnersTE
        
        SectionsProfiles = []
        SparLEFullWing=[]
        SparTEFullWing=[]
        SparMIDFullWing=[]


        SparLEup=[]
        SparLEdown=[]
        SparMIDup=[]
        SparMIDdown=[]
        SparTEup=[]
        SparTEdown=[]
        
        planeStiffEps=[]

  
        PStiffUpFullWing=[]
        PStiffDownFullWing=[]    
        StringerUpFullWing=[]
        StringerDownFullWing=[]

       
        EpsS= np.linspace(0, 1, self.SegmentNoLoft+1)                         
        SectionsProfiles = [self.AirfoilFunct(EpsS[d], LE[d], self.ChordFunct,
                                      ChordFactor, self.DihedralFunct,
                                      self.TwistFunct).Curve
                         for d in xrange(self.SegmentNoLoft+1)]    
        
# Stringers generation          
        h=0
        cont3=0
        PointsUp=[]
        PointsDown=[]
        PointsMIDup=[]
        PointsMIDdown=[]
        
        for SecR in SectionsProfiles:
            ChordLenght=((ChordFactor*self.ChordFunct(EpsS[h]))/
               np.cos(np.radians(self.TwistFunct(EpsS[h]))))
            LELengthX=LE[h,0]
            LEplanePositionX = (self.PositionLESpar*ChordLenght)*np.cos(self.TwistFunct(EpsS[h]))+LELengthX
            TEplanePositionX = (self.PositionTESpar*ChordLenght)*np.cos(self.TwistFunct(EpsS[h]))+LELengthX
            MIDplanePositionX = (self.PositionMIDSpar*ChordLenght)*np.cos(self.TwistFunct(EpsS[h]))+LELengthX
            
            planeLEi =Geom_Plane(gp_Pln(gp_Pnt(LEplanePositionX,0.,0.),gp_Dir(gp_Vec(gp_Pnt(0., 0., 0.),gp_Pnt(1.,0.,0.)))))
            planeMIDi =Geom_Plane(gp_Pln(gp_Pnt(MIDplanePositionX,0.,0.),gp_Dir(gp_Vec(gp_Pnt(0., 0., 0.),gp_Pnt(1.,0.,0.)))))
            planeTEi =Geom_Plane(gp_Pln(gp_Pnt(TEplanePositionX,0.,0.),gp_Dir(gp_Vec(gp_Pnt(0., 0., 0.),gp_Pnt(1.,0.,0.)))))

#            planeLEi.Rotate(gp_OY(), -np.radians(self.TwistFunct(EpsS[h])))
#            planeMIDi.Rotate(gp_OY(), -np.radians(self.TwistFunct(EpsS[h])))
#            planeTEi.Rotate(gp_OY(), -np.radians(self.TwistFunct(EpsS[h])))
            
            SegmentStiffLE=(MIDplanePositionX-LEplanePositionX)/(self.NoStiffnersLE-1)
            SegmentStiffTE=(TEplanePositionX-MIDplanePositionX)/(self.NoStiffnersTE)
            NoStiffners=self.NoStiffnersLE+self.NoStiffnersTE

            PStiffUpi=[]
            PStiffDowni=[]
            PosStiffiEps=LEplanePositionX
            
            for g in xrange(NoStiffners):
                
                planeStiffEps=Geom_Plane(gp_Pln(gp_Pnt(PosStiffiEps,0.,0.),gp_Dir(gp_Vec(gp_Pnt(0., 0., 0.),gp_Pnt(1.,0.,0.)))))
               
                StiffiEps = GeomAPI_IntCS(SecR, planeStiffEps.GetHandle())
                with act.assert_isdone(StiffiEps, 'failed to calculate intersection'):
                  nb_results = StiffiEps.NbPoints()
                  if nb_results == 1:
                    return StiffiEps.Point(1)
                   #print("One point intersection")
                  elif nb_results >= 1:
                   #print("More than one intersection point") 
                    PStiff = []
                    for i in range(1, nb_results+1):
                       PStiff.append(StiffiEps.Point(i))
                     
                  else:
                    return None
                PStiffUpi.append(PStiff[0])
                PStiffDowni.append(PStiff[nb_results-1])
                if g <=self.NoStiffnersLE-2:
                    PosStiffiEps=PosStiffiEps+SegmentStiffLE
                else:
                    PosStiffiEps=PosStiffiEps+SegmentStiffTE
                    
                PointsUp.append(PStiff[0])
                PointsDown.append(PStiff[nb_results-1])
                
                  
            SparLEi = GeomAPI_IntCS(SecR, planeLEi.GetHandle())
            with act.assert_isdone(SparLEi, 'failed to calculate intersection'):
               nb_results = SparLEi.NbPoints()
              #print(nb_results)
               if nb_results == 1:
                 return SparLEi.Point(1)
                #print("One point intersection")
               elif nb_results >= 1:
                #print("More than one intersection point")  
                 PLE = []
                 for i in range(1, nb_results+1):
                    PLE.append(SparLEi.Point(i))
                    
               else:
                 return None

            SegLE = GC_MakeSegment(PLE[0],PLE[nb_results-1]).Value() 
            
            SparMIDi = GeomAPI_IntCS(SecR, planeMIDi.GetHandle())
            with act.assert_isdone(SparMIDi, 'failed to calculate intersection'):
               nb_results = SparMIDi.NbPoints()
              #print(nb_results)
               if nb_results == 1:
                 return SparMIDi.Point(1)
                #print("One point intersection")
               elif nb_results >= 1:
                #print("More than one intersection point")  
                 PMID = []
                 for i in range(1, nb_results+1):
                    PMID.append(SparMIDi.Point(i))
                    
               else:
                 return None

            SegMID = GC_MakeSegment(PMID[0],PMID[nb_results-1]).Value()  
            PointsMIDup.append(PMID[0])
            PointsMIDdown.append(PMID[nb_results-1])
            
            SparTEi = GeomAPI_IntCS(SecR, planeTEi.GetHandle())
            with act.assert_isdone(SparTEi, 'failed to calculate intersection'):
               nb_results = SparTEi.NbPoints()
              #print(nb_results)
               if nb_results == 1:
                 return SparTEi.Point(1)
                #print("One point intersection")
               elif nb_results >= 1:
                #print("More than one intersection point")  
                 PTE = []
                 for i in range(1, nb_results+1):
                    PTE.append(SparTEi.Point(i))
                    
               else:
                 return None

            SegTE = GC_MakeSegment(PTE[0],PTE[nb_results-1]).Value()            

           #print(len(SparMIDdown))
           #print(len(SparMIDup))
           #print(len(SparLEup))
           #print(len(SparLEdown))
            SparLEFullWing.append(SegLE)
            SparMIDFullWing.append(SegMID)
            SparTEFullWing.append(SegTE)
            SparLEup.append(PLE[0])
            SparLEdown.append(PLE[nb_results-1])
            SparMIDup.append(PMID[0])
           #print(PMID[nb_results-1])
            SparMIDdown.append(PMID[nb_results-1])
            SparTEup.append(PTE[0])
            SparTEdown.append(PTE[nb_results-1])
            PStiffUpFullWing.append(PStiffUpi)
            PStiffDownFullWing.append(PStiffDowni)
            cont3=cont3+1
              
            h=h+1

# Method to create spars - VIA SEGMENTS
#LE     
        LoftLESparFullWing=[]                 
        for iii in xrange(0,len(PointsUp)-2*(self.NoStiffners-1),self.NoStiffners):
           SegUp=GC_MakeSegment(PointsUp[iii],PointsUp[iii+self.NoStiffners]).Value()    
           SegDown=GC_MakeSegment(PointsDown[iii+self.NoStiffners],PointsDown[iii]).Value() 
           SegSide1=GC_MakeSegment(PointsDown[iii],PointsUp[iii]).Value()    
           SegSide2=GC_MakeSegment(PointsUp[iii+self.NoStiffners],PointsDown[iii+self.NoStiffners]).Value()
           EdgeSegUp=act.make_edge(SegUp)
           EdgeSegDown=act.make_edge(SegDown)               
           EdgeSide1=act.make_edge(SegSide1)
           EdgeSide2=act.make_edge(SegSide2)
           WireSparLE=act.make_wire([EdgeSide1,EdgeSegUp,EdgeSide2,EdgeSegDown])
           FaceSparLEi=act.make_face(WireSparLE) 

           if self.ScaleFactor != 1:
                Origin = gp_Pnt(0., 0., 0.)
                FaceSparLE=act.scale_uniformal(FaceSparLEi, Origin, self.ScaleFactor)
                LoftLESparFullWing.append(FaceSparLE)   
#MID 
        LoftMIDSparFullWing=[]                 
        for iii in xrange(0,self.SegmentNoLoft):
           SegUp=GC_MakeSegment(PointsMIDup[iii],PointsMIDup[iii+1]).Value()    
           SegDown=GC_MakeSegment(PointsMIDdown[iii+1],PointsMIDdown[iii]).Value() 
           SegSide1=GC_MakeSegment(PointsMIDdown[iii],PointsMIDup[iii]).Value()    
           SegSide2=GC_MakeSegment(PointsMIDup[iii+1],PointsMIDdown[iii+1]).Value()
           EdgeSegUp=act.make_edge(SegUp)
           EdgeSegDown=act.make_edge(SegDown)               
           EdgeSide1=act.make_edge(SegSide1)
           EdgeSide2=act.make_edge(SegSide2)
           WireSparMID=act.make_wire([EdgeSide1,EdgeSegUp,EdgeSide2,EdgeSegDown])
           FaceSparMIDi=act.make_face(WireSparMID) 

           if self.ScaleFactor != 1:
                Origin = gp_Pnt(0., 0., 0.)
                FaceSparMID=act.scale_uniformal(FaceSparMIDi, Origin, self.ScaleFactor)
                LoftMIDSparFullWing.append(FaceSparMID)      

#TE
        LoftTESparFullWing=[]                 
        for iii in xrange(self.NoStiffners-1,len(PointsUp)-(self.NoStiffners),self.NoStiffners):
           SegUp=GC_MakeSegment(PointsUp[iii],PointsUp[iii+self.NoStiffners]).Value()    
           SegDown=GC_MakeSegment(PointsDown[iii+self.NoStiffners],PointsDown[iii]).Value() 
           SegSide1=GC_MakeSegment(PointsDown[iii],PointsUp[iii]).Value()    
           SegSide2=GC_MakeSegment(PointsUp[iii+self.NoStiffners],PointsDown[iii+self.NoStiffners]).Value()
           EdgeSegUp=act.make_edge(SegUp)
           EdgeSegDown=act.make_edge(SegDown)               
           EdgeSide1=act.make_edge(SegSide1)
           EdgeSide2=act.make_edge(SegSide2)
           WireSparTE=act.make_wire([EdgeSide1,EdgeSegUp,EdgeSide2,EdgeSegDown])
           FaceSparTEi=act.make_face(WireSparTE) 

           if self.ScaleFactor != 1:
                Origin = gp_Pnt(0., 0., 0.)
                FaceSparTE=act.scale_uniformal(FaceSparTEi, Origin, self.ScaleFactor)
                LoftTESparFullWing.append(FaceSparTE)                       
                    
        for g in range(self.NoStiffners):
           for gg in range(self.SegmentNoLoft):
              SegStringerUp=GC_MakeSegment(PointsUp[gg*self.NoStiffners+g],PointsUp[gg*self.NoStiffners+g+self.NoStiffners]).Value() 
              SegStringerDown=GC_MakeSegment(PointsDown[gg*self.NoStiffners+g],PointsDown[gg*self.NoStiffners+g+self.NoStiffners]).Value() 
              EdgeUp=act.make_edge(SegStringerUp)
              EdgeDown=act.make_edge(SegStringerDown)
              WireStringUp=act.make_wire(EdgeUp)
              WireStringDown=act.make_wire(EdgeDown)                      
              if self.ScaleFactor != 1:
                Origin = gp_Pnt(0., 0., 0.)
                StrUpFullWingi=act.scale_uniformal(WireStringUp, Origin, self.ScaleFactor)
                StrDownFullWingi=act.scale_uniformal(WireStringDown, Origin, self.ScaleFactor) 
               #print('hello')
             #print(self.ScaleFactor)
              StringerUpFullWing.append(StrUpFullWingi)
              StringerDownFullWing.append(StrDownFullWingi)    
        
      
              
        return LoftLESparFullWing, LoftTESparFullWing, LoftMIDSparFullWing, StringerUpFullWing, StringerDownFullWing, PointsUp, PointsDown
            
    def _BuildRibs(self, ChordFactor, ScaleFactor):
        """Generates Ribs surface along the wingspan, given the general,
        nondimensional parameters of the object (variations of chord length,
        dihedral, etc.) and the two scaling factors.

        Parameters
        ----------
        ChordFactor - float
            Factor againt which the standard shape is scaled in the chordwise
            direction
        ScaleFactor - float
            Factor againt which the standard shape is scaled in the world
            coordinates

        Returns
        -------
        Ribs : TopoDS_Shape
            The generated Ribs surface
        """
        SectionsRibs = []
        RibsFace = []
        

        Step=(self.SegmentNoLoft)/(self.NoRibs-1)
        LEPointsR = self._GenerateLeadingEdge()
        Eps=np.linspace(0,1,self.SegmentNoLoft+1)
        SectionsRibs = [self.AirfoilFunct(Eps[i], LEPointsR[i], self.ChordFunct,
                                  ChordFactor, self.DihedralFunct, self.TwistFunct).Curve
                   for i in xrange(0,self.SegmentNoLoft+1,Step)]             
        
        LESpar, TESpar, MIDSpar, StringerUp, StringerDown, PointsUp, PointsDown=self._BuildSpars(self.ChordFactor,self.ScaleFactor)
           
# RIBS generated from segments            
        for i in range(0,self.NoRibs):
            for ii in range(0,self.NoStiffners-1):
               SegUp=GC_MakeSegment(PointsUp[ii+i*self.NoStiffners*(Step)],PointsUp[ii+i*self.NoStiffners*(Step)+1]).Value()    
               SegDown=GC_MakeSegment(PointsDown[ii+i*self.NoStiffners*(Step)+1],PointsDown[ii+i*self.NoStiffners*(Step)]).Value() 
               SegSide1=GC_MakeSegment(PointsDown[ii+i*self.NoStiffners*(Step)],PointsUp[ii+i*self.NoStiffners*(Step)]).Value()    
               SegSide2=GC_MakeSegment(PointsUp[ii+i*self.NoStiffners*(Step)+1],PointsDown[ii+i*self.NoStiffners*(Step)+1]).Value()
               EdgeSegUp=act.make_edge(SegUp)
               EdgeSegDown=act.make_edge(SegDown)               
               EdgeSide1=act.make_edge(SegSide1)
               EdgeSide2=act.make_edge(SegSide2)
               WireRib=act.make_wire([EdgeSide1,EdgeSegUp,EdgeSide2,EdgeSegDown])
               FaceRibi=act.make_face(WireRib) 
               if self.ScaleFactor != 1:
                    Origin = gp_Pnt(0., 0., 0.)
                    FaceRib=act.scale_uniformal(FaceRibi, Origin, self.ScaleFactor)
                    RibsFace.append(FaceRib) 
                     
                
        return SectionsRibs, RibsFace  
        


    def _BuildPanels(self, ChordFactor, ScaleFactor):
        
        LESpar, TESpar, MIDSpar, StringerUp, StringerDown, PointsUp, PointsDown=self._BuildSpars(self.ChordFactor,self.ScaleFactor)
        PanelUp=[]
        PanelDown=[]


        cont=1
        for ii in range(0,len(PointsUp)-self.NoStiffners-1):
            if cont==self.NoStiffners:
                cont=0
            else:
                cont=cont
        
# Panel from segments
                SegH1=GC_MakeSegment(PointsUp[ii],PointsUp[ii+self.NoStiffners]).Value()
                SegH2=GC_MakeSegment(PointsUp[ii+1],PointsUp[ii+self.NoStiffners+1]).Value()
                SegV1=GC_MakeSegment(PointsUp[ii],PointsUp[ii+1]).Value()
                SegV2=GC_MakeSegment(PointsUp[ii+self.NoStiffners],PointsUp[ii+self.NoStiffners+1]).Value()
                EdgeHUp1=act.make_edge(SegH1)
                EdgeHUp2=act.make_edge(SegH2)               
                EdgeVUp1=act.make_edge(SegV1)
                EdgeVUp2=act.make_edge(SegV2)
                WirePanelUp=act.make_wire([EdgeHUp1,EdgeVUp1,EdgeHUp2,EdgeVUp2])
                PanelUpi=act.make_face(WirePanelUp)
# In this way the pannel will be ordered in alternance between upper and lower surface                
                SegH1=GC_MakeSegment(PointsDown[ii],PointsDown[ii+self.NoStiffners]).Value()
                SegH2=GC_MakeSegment(PointsDown[ii+1],PointsDown[ii+self.NoStiffners+1]).Value()
                SegV1=GC_MakeSegment(PointsDown[ii],PointsDown[ii+1]).Value()
                SegV2=GC_MakeSegment(PointsDown[ii+self.NoStiffners],PointsDown[ii+self.NoStiffners+1]).Value()
                EdgeHUp1=act.make_edge(SegH1)
                EdgeHUp2=act.make_edge(SegH2)               
                EdgeVUp1=act.make_edge(SegV1)
                EdgeVUp2=act.make_edge(SegV2)
                WirePanelDown=act.make_wire([EdgeHUp1,EdgeVUp1,EdgeHUp2,EdgeVUp2])
                PanelDowni=act.make_face(WirePanelDown)                
            
            
            
                if self.ScaleFactor!= 1:
                  Origin = gp_Pnt(0., 0., 0.)
                  PanelUpi= act.scale_uniformal(PanelUpi, Origin, self.ScaleFactor)
                  PanelDowni= act.scale_uniformal(PanelDowni, Origin, self.ScaleFactor)

            
                PanelUp.append(PanelUpi)
                PanelDown.append(PanelDowni)
            cont=cont+1

        return PanelUp, PanelDown     

        
    def GenerateStruct(self, ChordFactor, ScaleFactor):
        from OCC.TopoDS import TopoDS_Builder, TopoDS_Compound, TopoDS_Shape, TopoDS_HShape

        x0 = [ChordFactor, ScaleFactor]
                 
        SectionRibs,RibFace = self._BuildRibs(*x0)
        LSparLEFullWing, LSparTEFullWing, LSparMIDFullWing, StringerUpFullWing, StringerDownFullWing, PointsFullwingUp, PointsFullWingDown= self._BuildSpars(*x0)
        PanelUp, PanelDown=self._BuildPanels(*x0)
         
        builder=TopoDS_Builder()

# Add Ribs
        Ribs=TopoDS_Compound()   
        builder.MakeCompound(Ribs) 
        for g in RibFace:
           builder.Add(Ribs, g)
        self.AddComponent(Ribs,'Ribs')   
  
##Add Spars
        Spars=TopoDS_Compound()
        builder.MakeCompound(Spars)
        SparsLEMIDTE=LSparLEFullWing+LSparMIDFullWing+LSparTEFullWing
        for h in SparsLEMIDTE:    
           builder.Add(Spars, h)
        self.AddComponent(Spars,'Spars')     
#Add Stringers
        Stringers=TopoDS_Compound()    
        builder.MakeCompound(Stringers) 
        for p in xrange(self.NoStiffners*self.SegmentNoLoft):
          p1=StringerUpFullWing[p]
          p2=StringerDownFullWing[p]
          builder.Add(Stringers, p1)
          builder.Add(Stringers, p2)
        self.AddComponent(Stringers,'Stringers')      

#Add Panels
        Panels=TopoDS_Compound()
        builder.MakeCompound(Panels)
        for h in xrange(len(PanelUp)):
            h1=PanelUp[h]
            h2=PanelDown[h]
            builder.Add(Panels, h1)
            builder.Add(Panels, h2)
        self.AddComponent(Panels,'Panels')               
        vec = gp_Vec(gp_Pnt(-0.5, 0., 0.), self.ApexPoint)
        self.TranslateComponents(vec)           
        return None                                                 