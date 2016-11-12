#! /usr/bin/python
#Filename:   Autopilot_Dynamics_C_module.py
#Description:    Calculation module for Applications
#                                ____ Beta version ____

from MpltLib_GtkAggPlot import DrawDiag_GtkAgg as DrawDiag
import Tkinter as Tkn
import numpy as Np
from scipy.integrate import odeint
from scipy.signal import lti, impulse, step
from sys import exit

print 'Calculation module imported .........'

global const
global Omega
const=32.1768
Omega= Np.arange (1e-3, 1e3, 1e-3)


def DE_longitudinal (x, t):
    """ Description: Integrating function for longitudinal autopilot """
    return  [ Np.dot(x[0],A[0])+Np.dot(x[1],A[1])-Np.dot(x[2],const)+Np.dot(A[2],A[12]), 
                                                                                     Np.dot(x[0],A[3])+Np.dot(x[1],A[4])+Np.dot(x[3],A[11])+Np.dot(A[5],A[12]), 
                                                                                     x[3], 
                                                                                     Np.dot(x[0],A[6])+Np.dot(x[1],A[8])+Np.dot(A[10],
                                                                                                                                A[12])+Np.dot(x[3],A[9])+Np.dot(A[7],Np.dot(x[0],A[3])+Np.dot(x[1],
                                                                                                                                                                                              A[4])+Np.dot(x[3],
                                                                                                                                                                                              A[11])+Np.dot(A[5],A[12])) ]

def DE_lateral (x, t):
    """ Description: Integrating function for lateral autopilot """
    return [    Np.dot(x[0], B[0])+Np.dot(x[3], const/B[10])-x[2]+Np.dot(B[1], B[11]), 
                                                                              Np.dot(x[0], B[2])+Np.dot(x[1], B[3])+Np.dot(x[2], B[4])+Np.dot(B[5], B[11]),
                                                                              Np.dot(x[0], B[6])+Np.dot(x[1], B[7])+Np.dot(x[2], B[8])+Np.dot(B[9], B[11]), 
                                                                              x[1]  ]

def Autoplt_longitudinal (Longitudinal):
    print 'Longitudinal Autopilot'
    global A
    t0=0
    ImportData_Vars_long = ('tf','Xu','Xw','Xg','Zu','Zw','Zg','Mu','Mv','Mw','Mq','Mg','Uo','g')
    for i, value in enumerate (Longitudinal):
        if i==0:
            ImportData_Vals_long= []
        ImportData_Vals_long.append (value.get())
        
    for var, val in zip(ImportData_Vars_long, ImportData_Vals_long):
        print '%s = %5.3s' % (var, val)
        exec '%s=%s' % (var, val)
    print 'Test variables longitudinal      g=', g

    t= Np.arange (t0, tf, .01)
    x0= Np.zeros(4)
    A= Np.array([ Xu,Xw,Xg,Zu,Zw,Zg,Mu,Mv,Mw,Mq,Mg,Uo,g ])

    print 'A %s= %s' % (A.shape, A)
    print 't%s=%s' % (t.shape, t)
    x= odeint (DE_longitudinal, x0, t)
    print 'x%s=%s' % (x.shape, x)
    
    DrawDiag ([t], [x], FigureNum=1, NumSubplot=1, FigTitle='Longitudinal Autopilot', Xlabel=['Time (sec)'], Ylabel=[r'u, w, $\theta $'], render=True, geometry=(7, 6))

def Autoplt_lateral (Lateral):
    print 'Lateral Autopilot'
    global B
    t0=0
    ImportData_Vars_lat = ('tf','Yv', 'Yg_izv', 'Lbt', 'Lpt', 'Lrt', 'Lg_it', 'Nbt', 'Npt', 'Nvt', 'Ng_it', 'Uo', 'gi', 'g')
    for i, value in enumerate (Lateral):
        if i==0:
            ImportData_Vals_lat=[]
        ImportData_Vals_lat.append (value.get())
    
    for var, val in zip(ImportData_Vars_lat, ImportData_Vals_lat):
        print '%s = %5.3s' % (var, val)
        exec '%s=%s' % (var, val)
        
    print 'Test variables lateral       Lbt=', Lbt
    t= Np.arange(t0, tf, .01)
    x0= Np.zeros(4)
    B= Np.array ([ Yv, Yg_izv, Lbt, Lpt, Lrt, Lg_it, Nbt, Npt, Nvt, Ng_it, Uo, gi ])
    x= odeint (DE_lateral, x0, t)

    DrawDiag ([t], [x[:, :3]], FigureNum=1, NumSubplot=1, FigTitle='Lateral Autopilot', Xlabel=['Time (sec)'], Ylabel=[r'u, w, $\theta $'], render=True, geometry=(7, 6))
    
def Autoplt_transversal (Transversal):
    print 'Transversal Autopilot'
    ImportData_Vars_trans = ('tf', 'Yv', 'Yg_izv', 'Lbt', 'Lpt', 'Lrt', 'Lg_it', 'Nbt', 'Npt', 'Nvt', 'Ng_it', 'Uo', 'gi')
    
    for i, value in enumerate (Transversal):
        if i==0:
            ImportData_Vals_trans=[]
        ImportData_Vals_trans.append (value.get())
    
    for var, val in zip(ImportData_Vars_trans, ImportData_Vals_trans):
        print '%s = %5.3s' % (var, val)
        exec '%s=%s' % (var, val)
    print 'Test variables transversal       Lrt=', Lrt
    
    B= Np.array ([ Yv, Yg_izv, Lbt, Lpt, Lrt, Lg_it, Nbt, Npt, Nvt, Ng_it, Uo, gi ])
    A0= Np.array ([ [B[0], 0, -1, const/B[10]], 
                                           [B[2], B[3], B[4], 0],
                                           [B[6], B[7], B[8], 0], 
                                           [0, 1, 0, 0] ])
    B0= Np.array ([ [B[1]], [B[5]], [B[9]], [0] ])
    C0= Np.array ([ 1, 0, 0, 0]).reshape (1, 4)
    D0= Np.array ([0])
    
    GraphAgrams ([], [-1./227]*2, 1./227**2, A0, B0, C0, D0)
    
class GraphAgrams:
    def __init__ (self, Arg1, Arg2, Arg3, A, B, C, D):
        
        def Tf (*args, **kwargs):
            """ Description:
            Function for arranging  input array for Bode function plot  
            """
            if len(args) == 3:
                zeros= Np.asarray (args [0])
                poles= Np.asarray (args [1])
                gain= args [2]
                G= lti (zeros, poles, gain)
            elif len (args) == 2:
                B, A= args
                G= lti (B, A)
            G.desc= kwargs.get ('desc', ' ')
            print 'Tf function return parameter G=', G
            return G

        def Bode (G):
            '''     Description:    
            Plots Bode plot (Magnitude and Phase plot)   '''
            
            MpltLib.figure ()
            w= Np.arange (1e-4j, 1e-1j, 1e-6j)
            y= Np.polyval (G.num, w) / Np.polyval (G.den, w)
            mag= 20.*Np.log10(abs(y))
            pha= Np.arctan2 (Np.imag(y), Np.real(y)) *180. / Np.pi
            for i in Np.arange (1, pha.shape[0]):
                if abs (pha[i]-pha[i-1]) > 179:
                    pha[i:] -= 360.
            
            MpltLib.subplot (211)
            MpltLib.semilogx (w.imag, mag)
            MpltLib.grid ()
            MpltLib.ylabel (r'Magnitude (db)')
            
            MpltLib.subplot (212)
            MpltLib.semilogx  (w.imag, pha)
            MpltLib.grid ()
            MpltLib.ylabel (r'Phase(deg)')
            print 'Bode function return parameters mag=%s,pha=%s' % (mag, pha)        
            
            print 'pha =', pha.min
            print 'mag =', mag
            print "G=", G
            print "G.num=", G.num
            print "G.den=", G.den
            print "\t y=", y; print "\t y.imag=", Np.imag(y); print "\t y.real=", Np.real(y)
            print "w=", w; print 'w.shape, len(w) =', w.shape, len (w)
            print "w.imag=", w.imag; print Np.imag(w)
            print "A(%s)=%s" % (A.shape, A); print "B(%s)=%s" % (B.shape, B); print "C(%s)=%s" % (C.shape, C); print "D(%s)=%s" % (D.shape, D)
            
            return mag, pha
        
        def Impulse (G):
            MpltLib.figure ()
            MpltLib.plot (*G.impulse())
            MpltLib.grid()
            MpltLib.title ('Impulse responce: %s' %G.desc)
            MpltLib.xlabel ('Time (sec)')
            MpltLib.ylabel ('Amplitude')
        
        def Step (G):
            MpltLib.figure ()
            MpltLib.plot (*G.step())
            MpltLib.grid()
            MpltLib.title ('Step responce: %s ' %G.desc)
            MpltLib.xlabel ('Time (sec)')
            MpltLib.ylabel ('Amplitude')
        
        print  "A(%s)=%s\nB(%s)=%s\n,C(%s),=%s\nD(%s)=%s," % (A.shape, A, B.shape, B, C.shape, C, D.shape, D)
        Tff= lti (A, B, C, D)

        G = Tf ([], [-1./227]*2, 1./227**2)
        G.desc= r'$G_{H2}(s)=\frac{1}{(1+530s)^3)$'
        Bode (G)
        MpltLib.show ()
        Impulse (Tff)

        Step (Tff)
        MpltLib.show ()

class systems_info_v2:
    def __init__ (self, Arg1, Arg2, Arg3):
    
        def Tf (*args, **kwargs):
                '''     Description:
                Function for arranging  input array for Bode function plot  
                '''
                
                if len(args) == 3:
                    zeros= Np.asarray (args [0])
                    poles= Np.asarray (args [1])
                    gain= args [2]
                    G= lti (zeros, poles, gain)
                elif len (args) == 2:
                    B, A= args
                    G= lti (B, A)
                G.desc= kwargs.get ('desc', ' ')
                print 'Tf function return parameter G=', G
                return G
        
        def bodeplot (fi, f, tf, clear=True):
            MpltLib.figure (fi)
            if clear:
                MpltLib.clf ()
            
            MpltLib.subplot (211)
            MpltLib.semilogx (f, 20*Np.log10(abs(tf)))
            MpltLib.ylabel ("Mag.Ratio  (dB)")
            
            MpltLib.subplot (212)
            MpltLib.semilogx (f, arctan2(Np.imag(tf), Np.real(tf))*180./Np.pi)
            MpltLib.ylabel ("Phase (deg.)")
            MpltLib.xlabel ("Freg. (Hz)")

def LongTFFactors (InputData_Vars, InputData_Vals, Fig=None):
    print '\nLongitudinal Aerodynamics'    
    print InputData_Vars
    for i, value in enumerate (InputData_Vals):
        if i==0:
            InputData_ValZ=[]
        InputData_ValZ.append (value.get())
    for var, val in zip(InputData_Vars, InputData_ValZ):
        print '%s = %5.8s' % (var, val)
        exec '%s=%s' % (var, val)
    Theta_delta_num= Np.array ([
                            At,
                            At*(Tt1+Tt2), 
                            At*Tt1*Tt2  ])
    W_delta_num= Np.array ([
                        Aw, 
                        Aw*(2*Ksi_w*Omega_w+Tw1), 
                        Aw*(Omega_w**2+2*Ksi_w*Omega_w*Tw1), 
                        Aw*Omega_w**2*Tw1   ])
    H_delta_num= Np.array ([
                        Ah, 
                        Ah*(Th1+Th2+Th3), 
                        Ah*(Th1*Th2+(Th1+Th2)*Th3), 
                        Ah*Th1*Th2*Th3  ])
    Delta_long1= Np.array ([
                            1,
                            2*Tp1*Tp2,#Tp1+Tp2, 
                            Tp2**2 ])#Tp1*Tp2 ])
    Delta_long2= Np.array ([
                            1,
                            2*Tsp1*Tsp2,#Tsp1+Tsp2, 
                            Tsp2**2 ])#Tsp1*Tsp2   ])
    Delta_long= Np.polyval (Delta_long1, Omega)
    print 'Omega %s=%s' % (Omega.shape, Omega)
    print 'Poly1(%s)=%s' % (Np.polyval(Delta_long1, Omega).shape, Np.polyval(Delta_long1, Omega))
    print 'Delta_long (%s)=%s ' % (Delta_long.shape, Delta_long)
    
    Theta_delta_LTI= lti (Theta_delta_num, Delta_long)
    Theta_delta= Np.polyval (Theta_delta_LTI.num, Omega) / Np.polyval (Theta_delta_LTI.den, Omega)
    W_delta_LTI= lti (W_delta_num, Delta_long)
    W_delta= Np.polyval (W_delta_LTI.num, Omega) / Np.polyval (W_delta_LTI.den, Omega)
    H_delta_LTI= lti (H_delta_num, Delta_long)
    H_delta= Np.polyval (H_delta_LTI.num, Omega) / Np.polyval (H_delta_LTI.den, Omega)
    
    Theta_delta_Mag= 20.*Np.log10 (abs(Theta_delta))
    W_delta_Mag= 20.*Np.log10 (abs(W_delta))
    H_delta_Mag= 20.*Np.log10 (abs(H_delta))
    
    Theta_delta_Pha= Np.arctan2 (Np.imag(Theta_delta), Np.real(Theta_delta)) *180. / Np.pi
    W_delta_Pha= Np.arctan2 (Np.imag(W_delta), Np.real(W_delta)) *180. / Np.pi
    H_delta_Pha= Np.arctan2 (Np.imag(H_delta), Np.real(H_delta)) *180. / Np.pi
    
    MpltLib.figure (figsize=(7, 8), facecolor='w', edgecolor='w')              
    MpltLib.clf()
    MpltLib.suptitle ('Longitudinal Dynamics '+type, fontsize=12)
    MpltLib.suptitle ('\n\nBode Plot', fontsize=9)
    DrawDiag ([Omega]*2, [Theta_delta_Mag, Theta_delta_Pha], FigureNum=1, NumSubplot=2, FigTitle='Longitudinal Dynamics '+type,
               FigSuptitle='Bode Plot', Xlabel=['$Frequency (\omega (rad/sec))$']*2, Ylabel=['Magnitude (dB)', 'Phase (deg)'], type='semilogx',
               geometry=(7, 8))    
    
    MpltLib.figure (figsize=(7, 6), facecolor='w', edgecolor='w')                         
    MpltLib.clf()
    MpltLib.suptitle ('Longitudinal Dynamics '+type, fontsize=12)
    MpltLib.suptitle ('\n\nStep Plot', fontsize=9)
    DrawDiag ([Omega], [Theta_delta_LTI], FigureNum=2, NumSubplot=1, FigTitle='Longitudinal Dynamics '+type, FigSuptitle='Step Plot',
               Xlabel=['Time (sec)'], Ylabel=['Amplitude'], type='step')
    
    MpltLib.show ()

def LatTFFactors (InputData_Vars, InputData_Vals, Fig=None, type='unknown'):
    print '\nLateral Aerodynamics'+type
    print InputData_Vars
    for i, value in enumerate (InputData_Vals):
        if i==0:
            InputData_ValZ=[]
        InputData_ValZ.append (value.get())
    for var, val in zip(InputData_Vars, InputData_ValZ):
        print '%s = %5.8s' % (var, val)
        exec '%s=%s' % (var, val)
    
    A=1-Ixz**2/Ix/Iz
    
    Beta_delta_num= Np.array ([
                               Ab, 
                               Ab*(Tb1+Tb2+Tb3), 
                               Ab*(Tb1*Tb2+Tb3*(Tb1+Tb2)), 
                               Ab*Tb1*Tb2*Tb3   ])
    Delta_lat= Np.array ([
                               A, 
                               A*(2*Ksi_d*Omega_d+Ts+TR), 
                               A*((Omega_d**2+2*Ksi_d*Omega_d*(Ts+TR))+Ts*TR), 
                               A*(Omega_d**2*(Ts+TR)+2*Ksi_d*Omega_d*Ts*TR), 
                               A*(Omega_d**2*Ts*TR) ])
    Fi_delta_num= Np.array ([
                             Ap, 
                             Ap*2*Tp2*Tp3, 
                             Ap*Tp3 ])
    R_delta_num= Np.array ([
                            Ar, 
                            Ar*(2*Ksi_r*Omega_r+Tr1), 
                            Ar*(Omega_r**2+2*Ksi_r*Omega_r*Tr1), 
                            Ar*Omega_r**2*Tr1   ])
                
    Beta_delta_LTI= lti (Beta_delta_num, Delta_lat)
    Beta_delta= Np.polyval (Beta_delta_LTI.num, Omega) / Np.polyval (Beta_delta_LTI.den, Omega)
    Fi_delta_LTI= lti (Fi_delta_num, Delta_lat)
    Fi_delta= Np.polyval (Fi_delta_LTI.num, Omega)/Np.polyval (Fi_delta_LTI.den, Omega)
    R_delta_LTI= lti (R_delta_num, Delta_lat)
    R_delta= Np.polyval (R_delta_LTI.num, Omega)/Np.polyval (R_delta_LTI.den, Omega)
    
    print 'Beta_delta_LTI= ', Beta_delta_LTI.num, Beta_delta_LTI.den
    print 'Beta_delta= ', Beta_delta_num, Delta_lat

    Beta_delta_Mag= 20.*Np.log10 (abs(Beta_delta))
    Fi_delta_Mag= 20.*Np.log10 (abs(Fi_delta))
    R_delta_Mag= 20.*Np.log10 (abs(R_delta))
    
    Beta_delta_Pha= Np.arctan2 (Np.imag(Beta_delta), Np.real(Beta_delta)) *180. / Np.pi
    Fi_delta_Pha= Np.arctan2 (Np.imag(Fi_delta), Np.real(Fi_delta)) *180. / Np.pi
    R_delta_Pha= Np.arctan2 (Np.imag(R_delta), Np.real(R_delta)) *180. / Np.pi
    
    print 'Beta_delta_Mag= ', Beta_delta_Mag
    print 'Beta_delta_Pha= ', Beta_delta_Pha

    DrawDiag ([Omega]*2, [Beta_delta_Mag, Beta_delta_Pha], FigureNum=1, NumSubplot=2, FigTitle='Lateral Dynamics '+type,
               FigSuptitle='Bode Plot', Xlabel=['$Frequency (\omega (rad/sec))$']*2, Ylabel=['Magnitude (dB)', 'Phase (deg)'],
               type='semilogx', geometry=(7, 8), render=False)

    DrawDiag ([Omega], [Beta_delta_LTI], FigureNum=2, NumSubplot=1, FigTitle='Lateral Dynamics '+type, FigSuptitle='Step Plot',
               Xlabel=['Time (sec)'], Ylabel=['Amplitude'], type='step', render=True)
