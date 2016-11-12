#! /usr/bin/python
# FLIGHT DYNAMICS AUTOPILOT
#   Filename: Autopilot_Dynamics.py
#       Modules: Autopilot_Dynamics_Defmodule,
#                            Autopilot_Dynamics_C_module, 
#                                     ____ Beta version ____

from Autopilot_Dynamics_C_module import Autoplt_longitudinal, Autoplt_lateral, Autoplt_transversal
from Autopilot_Dynamics_C_module import LongTFFactors, LatTFFactors
import Autopilot_Dynamics_Def_module as Def
import Tkinter as Tkn
import numpy as Np

root= Tkn.Tk ()
root.geometry ('+250+20')
root.title ('System and Control of Aircraft')
App= Def.Application(root)

Canvs = Tkn.Canvas (root,  relief=Tkn.SUNKEN, bd=2, height=90, bg='white')
Canvs.pack (padx=2, pady=2, fill=Tkn.X)
Img = Tkn.PhotoImage (file='masinac.gif')
Canvs.create_image (70, 5, image=Img, anchor=Tkn.N)
Canvs.create_text (280, 50, text=u''' University of Belgrade - Faculty of Mechanical Engineering
Chair of Aeronautics and Astronautics
Course title: Systems and Control of Aircraft''',  font=Def.Font9B, justify=Tkn.CENTER )

Lbl_title = Tkn.Label (root
                       , text=u'''Program for calculation and graphical presentation of Lateral and Longitudinal Dynamics 
- Autopilot Characteristics''', font=Def.Font9B)
Lbl_title.pack(padx=10, pady=5)


##Longitudinal Autopilot
Tlevel_longitudinal = Tkn.Toplevel (root)
Tlevel_longitudinal.title ('Longitudinal autopilot')
Tlevel_longitudinal.geometry ('+250+260')
ImportData_Vars_long = ('tf','Xu','Xw','Xg','Zu','Zw','Zg','Mu','Mv','Mw','Mq','Mg','Uo','g')
ImportData_Vals_long = (10., .0097, .0016, 0., .0955, -1.43, -69.8, 0, -.0013, -.0235, -1.92, -26.1, 660., 32.17)

Longitudinal = Def.Multiple_Entries ( Tlevel_longitudinal, zip(ImportData_Vars_long, ImportData_Vals_long)
                  , Np.zeros(len(ImportData_Vars_long))+10, 'Input data \nlongitudinal autopilot  '  )
CBttn1 = Tkn.Button (Tlevel_longitudinal, text='Calculate parameters', font=Def.Font8, 
                     command=(lambda:Autoplt_longitudinal(Longitudinal)))
CBttn1.pack(padx=6, pady=3)


##Lateral Autopilot
Tlevel_lateral = Tkn.Toplevel (root)
Tlevel_lateral.title ('Lateral autopilot')
Tlevel_lateral.geometry ('+410+260')
ImportData_Vars_lat = ('tf','Yv', 'Yg_izv', 'Lbt', 'Lpt', 'Lrt', 'Lg_it', 'Nbt', 'Npt', 'Nvt', 'Ng_it', 'Uo', 'gi')
ImportData_Vals_lat = (10., -.1008, 0, -2.71, -1.232, .397, -1.62, 1.301, .0346, .257, .01875, 468.2, 32.17)

Lateral = Def.Multiple_Entries (Tlevel_lateral, zip(ImportData_Vars_lat, ImportData_Vals_lat)
                  , Np.zeros(len(ImportData_Vars_lat))+10, 'Input data \nlateral autopilot  ' )
CBttn2 = Tkn.Button (Tlevel_lateral, text='Calculate parameters', font=Def.Font8, 
                     command=(lambda:Autoplt_lateral(Lateral)))
CBttn2.pack(padx=6, pady=3)


##Transversal Autopilot
#   New approach - Examples from McRuer book
##LateraL Autopilot (McRuer example1 Table6.1 coefficients, Lateral Dimensional Derivatives)
Tlevel_Lateral1= Def.AppWindow (root, 'Lateral Dynamics', '+10+40')
LateralGIparameters_EntVars= ('Height (ft)', 'Mach number', u'Air Dencity (slugs/ft\xb3)', 'Airspeed (ft/sec)', 'Weight (lb)', 'Mass (slugs)', u'Ix (slug-ft\xb2)', u'Iz (slug-ft\xb2)' ,u'Ixz (slug-ft\xb2)')
LateralGIparameters_Vars= ('h', 'M0', 'ro', 'U0', 'W', 'm', 'Ix', 'Iz', 'Ixz')
LateralGIparameters_EntVals= (15e3, .4, .001496, 423, 17578, 546, 8190, 29140, 1994 )

LateralDimensDerivatives_EntVars= (u'Yv (1/sec)', u'Y\u03B4\u2090 (ft/sec\xb2)/rad', u'Y\u03B4\u2090* (1/sec)/rad', u'Y\u03B4r (ft/sec\xb2)/rad', u'Y\u03B4r* (1/sec)/rad',
u'L\u03B2 (1/sec\xb2)', u'Lp (1/sec)', 'Lr (1/sec)', u'L\u03B4\u2090 (1/sec\xb2)', u'L\u03B4r (1/sec\xb2)', 
u'N\u03B2 (1/sec\xb2)', 'Np (1/sec)', 'Nr (1/sec)', u'N\u03B4\u2090 (1/sec\xb2)', u'N\u03B4r (1/sec\xb2)')
LateralDimensDerivatives_Vars= ('Yv', 'Yda','YdaT', 'Ydr', 'YdrT', 'Lb', 'Lp', 'Lr', 'Lda', 'Ldr', 'Nb', 'Np', 'Nr', 'Nda', 'Ndr')
LateralDimensDerivatives_EntVals= (.1476, .892, -.00211, 8.92, .0211, 14.02, .987, .608, 8.76, 2.8, 8.21, 0, .4, .1641, 3.64)

LongitudinalDimensDerivatives_EntVars= ('Xu', 'Xw', u'X\u03B4\u2091', 'Zu', 'Zw', u'Z\u03B4\u2091', 'Mu', 'Mw', 'Mw*', 'Mq', u'M\u03B4\u2091')
LongitudinalDimensDerivatives_Vars= ('Xu', 'Xw', 'Xde', 'Zu', 'Zw', 'Zde', 'Mu', 'Mw', 'MwT', 'Mq', 'Mde')
LongitudinalDimensDerivatives_EntVals= (-.0097, .0016, 0, -.0955, -1.43, -69.8, 0, -.0235, -0.0013, -1.92, -26.1)

LateraL1 = Def.Multiple_Entries  (Tlevel_Lateral1, zip(LateralGIparameters_EntVars, LateralGIparameters_EntVals)
                    , Np.zeros(len(LateralGIparameters_EntVars))+10, 'Geometrical and Inertial Parameters') 
LateraL2= Def.Multiple_Entries (Tlevel_Lateral1, zip(LateralDimensDerivatives_EntVars, LateralDimensDerivatives_EntVals), 
                                Np.zeros(len(LateralDimensDerivatives_EntVars))+10, 'Lateral Dimensional Derivatives')
LateraL3= Def.Multiple_Entries (Tlevel_Lateral1, zip(LongitudinalDimensDerivatives_EntVars, LongitudinalDimensDerivatives_EntVals), 
                                Np.zeros (len(LongitudinalDimensDerivatives_EntVars))+10, 'Longitudinal Dimensional Derivatives')


##Lateral Autopilot (McRuer example2 Lateral and Longitudinal Transfer Function Derivatives)
class LateralTFFactors_aileron (Tkn.Toplevel):    
    def __init__ (self, parent=None):
        Tkn.Toplevel.__init__(self, parent)
        Tkn.Toplevel.title (self,'Lateral Aerodynamics')
        Tkn.Toplevel.geometry (self,'+580+260')
        LblFrm= Tkn.LabelFrame (self, text='Aileron Lateral Transfer\nFunction Factors', font='Helvetica 8 bold',
                                       relief=Tkn.SUNKEN, bd=1, labelanchor=Tkn.NE)
        LblFrm.pack(side=Tkn.TOP, padx=4, pady=4, anchor=Tkn.W, fill=Tkn.X)
        self.DLat= (('1/Ts', '1/TR', u'\u03B6d', u'\u03C9d'), 
                    ('Ts', 'TR', 'Ksi_d', 'Omega_d'), 
                    (.00509, 1.601, .0949, 3.06))
        self.NdFi= (('Ap', '1/Tp1', u'\u03B6p', u'\u03C9p'), 
                    ('Ap', 'Tp1', 'Tp2', 'Tp3'), 
                    (8.95, 0, .0954, 2.82))
        self.Ndr= (('Ar', '1/Tr1', u'\u03B6r', u'\u03C9r'), 
                   ('Ar', 'Tr1', 'Ksi_r', 'Omega_r'), 
                   (.777, 1.778, .548, 1.975))
        self.Ndb= ((u'A\u03B2', u'1/T\u03B21', u'1/T\u03B22 (\u03B6b)', u'1/T\u03B23 (\u03C9b)'), 
                                                                                      ('Ab', 'Tb1', 'Tb2', 'Tb3'), 
                                                                                      (-.00211, 368, .935, .583))
        Lat1= Def.Multiple_Entries (LblFrm, zip(self.DLat[0], self.DLat[2]), Np.ones(4)+10, u'\u0394lat', 'Helvetica 10 bold')
        Lat2= Def.Multiple_Entries (LblFrm, zip(self.NdFi[0], self.NdFi[2]), Np.ones(4)+10, u'N\u03B4\u03C6') 
        Lat3= Def.Multiple_Entries (LblFrm, zip(self.Ndr[0], self.Ndr[2]), Np.ones(4)+10, u'N\u03B4r')
        Lat4= Def.Multiple_Entries (LblFrm, zip(self.Ndb[0], self.Ndb[2]), Np.ones(4)+10, u'N\u03B4b') 
        CBttn = Tkn.Button (self, text='Calculate parameters', font=Def.Font8, 
                            command=(lambda: LatTFFactors (self.DLat[1]+self.NdFi[1]+self.Ndr[1]+self.Ndb[1]+LateralGIparameters_Vars,
                                                            Lat1+Lat2+Lat3+Lat4+LateraL1, Fig=1, type='-aileron')))
        CBttn.pack(padx=6, pady=3)

LateralTFFactors_aileron()

class LateralTFFactors_rudder (Tkn.Toplevel):
    def __init__ (self, parent=None):
        Tkn.Toplevel.__init__ (self, parent)
        Tkn.Toplevel.title (self, 'Lateral Dynamics')
        Tkn.Toplevel.geometry (self, '+790+260')    
        LblFrm= Tkn.LabelFrame (self, text='Rudder Lateral Transfer\n Function Factors', font='Helvetica 8 bold',
                                       relief=Tkn.SUNKEN, bd=1, labelanchor=Tkn.NE)
        LblFrm.pack(side=Tkn.TOP, padx=4, pady=4, anchor=Tkn.W, fill=Tkn.X)                                   
        self.DLat= (('1/Ts', '1/TR', u'\u03B6d', u'\u03C9d'), 
                    ('Ts', 'TR', 'Ksi_d', 'Omega_d'), 
                    (.00509, 1.016, .0949, 3.06))
        self.NdFi= (('Ap', '1/Tp1', '1/Tp2', '1/Tp3'), 
                    ('Ap', 'Tp1', 'Tp2', 'Tp3'), 
                    (3.75, 0, -2.89, 2.65))
        self.Ndr= (('Ar', '1/Tr1', u'\u03B6r', u'\u03C9r'), 
                   ('Ar', 'Tr1', 'Ksi_r', 'Omega_r'), 
                   (-3.9, 1.295, .198, .656))
        self.Ndb= ((u'A\u03B2', u'1/T\u03B21', u'1/T\u03B22', u'1/T\u03B23'), 
                   ('Ab', 'Tb1', 'Tb2', 'Tb3'), 
                   (.0211, -.021, 1.032, 185.4))
        Lat1= Def.Multiple_Entries (LblFrm, zip(self.DLat[0], self.DLat[2]), Np.ones(4)+10, u'\u0394lat', 'Helvetica 10 bold')
        Lat2= Def.Multiple_Entries (LblFrm, zip(self.NdFi[0], self.NdFi[2]), Np.ones(4)+10, u'N\u03B4\u03C6') 
        Lat3= Def.Multiple_Entries (LblFrm, zip(self.Ndr[0], self.Ndr[2]), Np.ones(4)+10, u'N\u03B4r')
        Lat4= Def.Multiple_Entries (LblFrm, zip(self.Ndb[0], self.Ndb[2]), Np.ones(4)+10, u'N\u03B4b')
        CBttn = Tkn.Button (self, text='Calculate parameters', font=Def.Font8,
                            command=(lambda: LatTFFactors(self.DLat[1]+self.NdFi[1]+self.Ndr[1]+self.Ndb[1]+LateralGIparameters_Vars,
                                                           Lat1+Lat2+Lat3+Lat4+LateraL1, Fig=3, type='-rudder')))
        CBttn.pack(padx=6, pady=3)

LateralTFFactors_rudder()

class LongitudinalTFFactors_elevator (Tkn.Toplevel):
    def __init__  (self, parent=None):
        Tkn.Toplevel.__init__ (self, parent)
        Tkn.Toplevel.title (self, 'Longitudinal Dynamics')
        Tkn.Toplevel.geometry (self, '+970+260')
        LblFrm= Tkn.LabelFrame (self, text='Elevator Llongitudinal Transfer\n Function Factors', font='Helvetica 8 bold',
                                       relief=Tkn.SUNKEN, bd=1, labelanchor=Tkn.NE)
        LblFrm.pack(side=Tkn.TOP, padx=4, pady=4, anchor=Tkn.W, fill=Tkn.X)                                   
        self.DLong= ((u'\u03B6sp (1/Tsp1)', u'\u03C9sp (1/Tsp2)', u'\u03B6p (1/Tp1)', u'\u03C9p (1/Tp2)'), 
                                                                                                 ('Tsp1', 'Tsp2', 'Tp1', 'Tp2'), 
                                                                                                 (.286, 2.45, .0439, .098))
        self.NdTheta= ((u'A\u03B8', u'1/T\u03B81', u'1/T\u03B82'),
                       ('At', 'Tt1', 'Tt2'), 
                       (7.4, .00308, .49) )
        self.Ndw= (('Aw', '1/Tw1', u'\u03B6w', u'\u03C9w'), 
                   ('Aw', 'Tw1', 'Ksi_w', 'Omega_w'), 
                   (22.9, 137.1, .069, .1067))
        self.Ndh= (('Ah', '1/Th1', '1/Th2', '1/Th3'), 
                   ('Ah', 'Th1', 'Th2', 'Th3'), 
                   (22.9, -7.72, .0208, 8.44))
        Long1= Def.Multiple_Entries (LblFrm, zip(self.DLong[0], self.DLong[2]), Np.ones(4)+10, u'\u0394long', 'Helvetica 10 bold')
        Long2= Def.Multiple_Entries (LblFrm, zip(self.NdTheta[0], self.NdTheta[2]), Np.ones(3)+10, u'Nd\u03B8')
        Long3= Def.Multiple_Entries (LblFrm, zip(self.Ndw[0], self.Ndw[2]), Np.ones(4)+10, u'N\u03B4w')
        Long4= Def.Multiple_Entries (LblFrm, zip(self.Ndh[0], self.Ndh[2]), Np.ones(4)+10, u'N\u03B4h')
        CBttn = Tkn.Button (self, text='Calculate parameters', font=Def.Font8, 
                            command=(lambda:LongTFFactors(self.DLong[1]+self.NdTheta[1]+self.Ndw[1]+self.Ndh[1], Long1+Long2+Long3+Long4)))
        CBttn.pack(padx=6, pady=3)

root.mainloop ()
