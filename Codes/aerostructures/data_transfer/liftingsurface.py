# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 15:52:58 2015

 LIFTINGSURFACE.PY ============================================================
 This module contains the definition of the class of 3d lifting surfaces.
 This class can be instantiated to generate wings, tailplanes, fins, propeller-
 or rotor blades, etc.

 This is an OCC_AirCONICS file, based on the Rhino 'AirCONICS' plugin
 by A. Sobester: https://github.com/sobester/AirCONICS
 ==============================================================================

@author: pchambers
"""
import numpy as np
from airconics.base import AirconicsShape
import AirCONICStools as act

from OCC.gp import gp_Pnt, gp_Vec
from OCC.GeomAbs import GeomAbs_C2
from OCC.GeomAbs import GeomAbs_Line

class LiftingSurface(AirconicsShape):
    """Airconics class for defining lifting surface shapes

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

    AirfoilFunct - function
        function defining the sectional Airfoil (see primitives.Airfoil)
        vs epsilon spanwise variable coordinate between 0 and 1
        (curvilinear attached)

    ChordFactor - int (default = 1)
        Scaling factor applied in chordwise direction

    ScaleFactor - int (default = 1)
        Scaling factor applied in all directions (uniform)

    OptimizeChordScale - int or bool (default = 0)
        TODO: Not yet used.
        
    LooseSurf - (default = 1)
        TODO: 
        
    SegmentNo - int (default = 11)
        Number of segments to sample the wing defined by input functions
    
    TipRequired - bool (default = False)
        TODO: Not yet used
        adds the wing tip face to components if true
    
    max_degree - (default = 8)
        maximum degree of the fitted NURBS surface

    continuity - OCC.GeomAbs.GeomAbs_XX Type
        the order of continuity i.e. C^0, C^1, C^2... would be 
        GeomAbs_C0, GeomAbs_C1, GeomAbs_C2 ...
        
    Attributes
    ----------
    self['Surface'] : TopoDS_Shape
        The generated lifting surface

    Notes
    -----
    * It is expected that users will create shapes mostly on initialisation
      of a LiftingSurface instance. GenerateLiftingSurface is therefore not
      expected to be called directly.
    
    * Output surface is stored in self['Surface']

    * See airconics.examples.wing_example_transonic_airliner for 
      example input functions
    
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
                 SegmentNo=10,
                 NPaero_chord=1,
                 NPaero_span=1,
                 TipRequired=False,
                 max_degree=8,
                 continuity=GeomAbs_C2                 
                 ):
#        Initialise the components using base class:
        super(LiftingSurface, self).__init__(components={},
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
                                         SegmentNo=SegmentNo,
                                         NPaero_chord=NPaero_chord,
                                         NPaero_span=NPaero_span,
                                         TipRequired=TipRequired,
                                         _max_degree=max_degree,
                                         _Cont=continuity
                                         )




        self.GenerateLiftingSurface(self.ChordFactor, self.ScaleFactor)

    def _GenerateLeadingEdge(self):
        """Epsilon coordinate attached to leading edge defines sweep
         Returns airfoil leading edge points
         """
        SegmentLength = 1. / self.SegmentNo

#       Array of epsilon at segment midpoints (will evaluate curve here)
        Epsilon_midpoints = np.linspace(SegmentLength/2., 1-SegmentLength/2.,
                                        self.SegmentNo)

#       We are essentially reconstructing a curve from known slopes at
#       known curve length stations - a sort of Hermite interpolation
#       without knowing the ordinate values. If SegmentNo -> Inf, the
#       actual slope at each point -> the sweep angle specified by
#       SweepFunct
        Tilt_array = self.DihedralFunct(Epsilon_midpoints)
        Sweep_array = self.SweepFunct(Epsilon_midpoints)

        DeltaXs =  SegmentLength * np.sin(Sweep_array*(np.pi/180.))
        DeltaYs = (SegmentLength * np.cos(Tilt_array*np.pi/180.) *
                                   np.cos(Sweep_array*np.pi/180.))
        DeltaZs = DeltaYs*np.tan(Tilt_array*np.pi/180.)

#        Initialise LE coordinate arrays and add first OCC gp_pnt at [0,0,0]:
#        Note: Might be faster to bypass XLE arrays and use local x only
#        XLE = np.zeros(self.SegmentNo + 1)
#        YLE = np.zeros(self.SegmentNo + 1)
#        ZLE = np.zeros(self.SegmentNo + 1)
#        LEPoints = [gp_Pnt(XLE[0], YLE[0], ZLE[0])]
        LEPoints = np.zeros((self.SegmentNo + 1, 3))

#        for i in xrange(self.SegmentNo):
#            XLE[i+1] = XLE[i] + DeltaXs[i]
#            YLE[i+1] = YLE[i] + DeltaYs[i]
#            ZLE[i+1] = ZLE[i] + DeltaZs[i]
#            LEPoints[i+1, :] = XLE[i+1], YLE[i+1], ZLE[i+1]
        Deltas = np.vstack([DeltaXs, DeltaYs, DeltaZs]).T
        LEPoints[1:, :] = np.cumsum(Deltas, axis=0)

        return LEPoints

    def _BuildLS(self, ChordFactor, ScaleFactor):
        """Generates a tentative lifting surface, given the general,
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
        LS : TopoDS_Shape
            The generated Lifting surface

        ActualSemiSpan : scalar
            TODO: currently not calculated, None is returned

        LSP_area : scalar
            TODO: currently not calculated, None is returned

        AR : scalar
            TODO: currently not calculated, None is returned

        WingTip : TopoDS face or shape
            TODO: currently not calculated, None is returned
        """
        LEPoints = self._GenerateLeadingEdge()

        Sections = []
        # TODO: These lists are used for when the curve has been smoothed or
        # the loft has failed, neither of which have been implemented yet
#        ProjectedSections = []
#        TEPoints_u = []
#        TEPoints_l = []

        Eps = np.linspace(0, 1, self.SegmentNo+1)
        Sections = [self.AirfoilFunct(Eps[i], LEPoints[i], self.ChordFunct,
                                      ChordFactor, self.DihedralFunct,
                                      self.TwistFunct).Curve
                    for i in xrange(self.SegmentNo+1)]

        self._Sections = Sections   # I used from debugging - remove it?

        # TODO: Implement chord projection and Curve start/end points 
        # to rescale smoothed curves and for secondary loft methods

        LS = act.AddSurfaceLoft(Sections,
                                max_degree=self._max_degree,
                                continuity=self._Cont)
                   

        if LS is None:
            pass


        WingTip = None

        if self.TipRequired:
            pass

        # Scaling
        if self.ScaleFactor != 1:
            Origin = gp_Pnt(0., 0., 0.)
            LS = act.scale_uniformal(LS, Origin, self.ScaleFactor)
            # TODO: Wing tip scaling (TipRequired is not implemented yet)
            if self.TipRequired and WingTip:
                pass

        self.RootChord = (self.ChordFunct(0)*ChordFactor)*ScaleFactor

        # Temporarily set other variables as None until above TODO's are done
        ActualSemiSpan = None
        LSP_area = None
        AR = None

        return LS, ActualSemiSpan, LSP_area, AR, WingTip
        
    def _BuildSurfacePoints(self, ChordFactor, ScaleFactor):
         LEPoints = self._GenerateLeadingEdge()

         SurfacePntsUp=[]
         SurfacePntsDown=[]         
         Sections=[]
#         SurfacePntstoAdd=[]
         
         Eps = np.linspace(0, 1, self.NPaero_span+1)
         Sections = [self.AirfoilFunct(Eps[i], LEPoints[i], self.ChordFunct,
                                      ChordFactor, self.DihedralFunct,
                                      self.TwistFunct).Curve
                    for i in xrange(self.NPaero_span+1)]
         numb=0
        
         for num in xrange(self.NPaero_span+1):
            section=Sections[num]
            sectionPnt=act.Uniform_Points_on_Curve(section, self.NPaero_chord)
            numb=numb+1
            for numm in xrange(1,self.NPaero_chord+1):
                print(sectionPnt[numm-1].Y())
                print(sectionPnt[numm-1].X())
                edgePnts=act.make_vertex(sectionPnt[numm-1])
                if self.ScaleFactor!= 1:
                  Origin = gp_Pnt(0., 0., 0.)   
                  Pnts=act.scale_uniformal(edgePnts, Origin, self.ScaleFactor)
                  if numm <= (self.NPaero_chord-1)/2+1:
                      SurfacePntsUp.append(Pnts)
                      if numm==(self.NPaero_chord-1)/2+1:
                          SurfacePntsDown.append(Pnts)
                  elif numm > (self.NPaero_chord-1)/2+1:
#                      SurfacePntsDown.append(SurfacePntstoAdd)
                      SurfacePntsDown.append(Pnts)
#                      SurfacePntsDown=SurfacePntsD+SurfacePntstoAdd
#            numbb=numbb+1    
#         if self.ScaleFactor!= 1:
#             Origin = gp_Pnt(0., 0., 0.)   
#             SurfacePnts=act.scale_uniformal(SurfacePnts, Origin, self.ScaleFactor) 

         print(np.shape(SurfacePntsUp))    
         print(np.shape(SurfacePntsDown)) 
         return SurfacePntsUp, SurfacePntsDown



    def GenerateLiftingSurface(self, ChordFactor, ScaleFactor):
        from OCC.TopoDS import TopoDS_Builder, TopoDS_Compound, TopoDS_Shape, TopoDS_HShape
        """Builds a lifting surface (wing, tailplane, etc.) with the Chord and
        Scale factors defined in inputs
        
        If OptimizeChordScale was specified on construction of this 
        LiftingSurface class, an optimized ChordFactor and ScaleFactor is found
        instead, with the local search started from the two given values.
        
        Parameters
        ----------
        ChordFactor : scalar
            The scaling factor to apply in the chordwise direction
        
        ScaleFactor : scalar
            the scaling factor to apply uniformly in all directions
        
        Returns
        -------
        None
        
        Notes
        -----
        Called on initialisation of a lifting surface class. Adds a
        ('Surface': Shape) key value pair to self.
        
        :Example:
            >>> Wing = liftingsurface.LiftingSurface(P,
                                                mySweepAngleFunction, 
                                                myDihedralFunction, 
                                                myTwistFunction, 
                                                myChordFunction, 
                                                myAirfoilFunction)
            >>> Surface = Wing['Surface']
        
        See Also
        --------
        airconics.examples.wing_example_transonic_airliner
        """
        x0 = [ChordFactor, ScaleFactor]

        LS, ActualSemiSpan, LSP_area, AR, WingTip = \
            self._BuildLS(*x0)
        
        PointsSurfUp, PointsSurfDown=self._BuildSurfacePoints(*x0)
#        Update instance components:
        self.AddComponent(LS, 'SurfaceLoft')
        builder=TopoDS_Builder()
        PointsUp=TopoDS_Compound()   
        builder.MakeCompound(PointsUp) 
        PointsDown=TopoDS_Compound()   
        builder.MakeCompound(PointsDown) 
        for g in PointsSurfUp:
           builder.Add(PointsUp, g)
        self.AddComponent(PointsUp,'PointsSurfUp') 
        for gg in PointsSurfDown:
           builder.Add(PointsDown, gg)
        self.AddComponent(PointsDown,'PointsSurfDown') 


        # Position the Components at the apex:
        vec = gp_Vec(gp_Pnt(0., 0., 0.), self.ApexPoint)
        self.TranslateComponents(vec)
        
        return None
