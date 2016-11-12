#! /usr/bin/python
#Filename:   Autopilot_Dynamics_Def_module.py
#Description:   Definition module for application
#                       ____ Beta version ____

import Tkinter as Tkn
import sys
import numpy as Np

print 'Definition module imported ....'

def TWindow (self, Lev1, Tlt1, Geom1):            
    Tlevel = Tkn.Toplevel (Lev1)
    Tlevel.title (Tlt1)
    Tlevel.geometry (Geom1)        
    return Tlevel

class AppWindow (Tkn.Toplevel):
    def __init__ (self, parent=None, Tlt='No title', Geom='+160+260'):
        Tkn.Toplevel.__init__ (self, parent)
        Tkn.Toplevel.title (self, Tlt)
        Tkn.Toplevel.geometry (self, Geom)
        self.parent= parent

class Application (Tkn.Frame):
    ''' Description:
    This class is used to make main window and Quit Button'''
    def __init__ (self, parent=None):
        Tkn.Frame.__init__(self, parent)
        self.parent=parent
        self.pack(side=Tkn.BOTTOM, padx=6, pady=7, anchor=Tkn.S)
        
        QBttn=Tkn.Button (self, text='Quit Application', command=(lambda: self.Exit()), padx=4, pady=2, width=20,   
                                                                                   font='Helvetica 8 bold',state=Tkn.ACTIVE)        
        QBttn.pack(anchor=Tkn.S)
    
    def Exit (self):
        print 'Application terminated ......'
        sys.exit()
    
def Multiple_Entries (self, InputArray, EntWidth, Labeltitle='No title', font='Helvetica 8 bold'
                      , SideFormat=Tkn.TOP, Lformat=Tkn.TOP):
    ''' Description:
    Function for creating multiple Entry fields for application GUI
    '''
    LblFrm= Tkn.LabelFrame (self, text=Labeltitle, font=font, relief=Tkn.SUNKEN, bd=1,
                             labelanchor=Tkn.NE)
    LblFrm.pack(side=Lformat, padx=4, pady=4, anchor=Tkn.W, fill=Tkn.X)
    
    for i, Data in enumerate(InputArray):
        if i==0:        
            OutputArray=[]
        Frm=Tkn.Frame(LblFrm)
        Frm.pack(side=SideFormat, fill=Tkn.X)
        Lbl=Tkn.Label(Frm, text=Data[0], font='Helvetica 8 ', justify=Tkn.LEFT)
        Lbl.pack(side=Tkn.LEFT, anchor=Tkn.W, padx=5, fill=Tkn.X)
        EntVar=Tkn.StringVar()
        Ent=Tkn.Entry (Frm, textvariable=EntVar, width=int(EntWidth[i]), justify=Tkn.CENTER)
        Ent.pack(side=Tkn.RIGHT, anchor=Tkn.E, padx=5, fill=Tkn.X)
        EntVar.set(Data[1])
        OutputArray.append(EntVar)
    return OutputArray

'''Description: 
Font library'''
Font8B= 'Helvetica 8 bold'
Font8= 'Helvetica 8'
Font9B= 'Helvetica 9 bold'
Font9= 'Helvetica 9'


def DebugInfo (VarArray, ValArray, Checkpoint, Status='Unknown'):
    '''     Description:
    Function for debugging application
    '''
    print '-----------------\nDebug Info\n'
    print 'Checkpoint %s \n' % Checkpoint
    for element, value in zip (VarArray, ValArray):
        try:
            print '%s (%s) = %s' % (element, Np.shape(value), value)
        except NameError:
            print '%s element is missing' % element
    print 'Checkpoint %s >> %s' % (Checkpoint, Status)
    print '-----------------'

def OutputFormatData (ValArray, Linetitle, Formatpattern, adjust=[2, 2]):
    '''     Description:
    Function for printing formated application output
    '''
    for n, row in enumerate (ValArray):
        for num, column in enumerate(row):
            print  Formatpattern.format (Linetitle[n].ljust(adjust[0]), 
                                                                   str(column).ljust(adjust[1])), 
            if num==len(row)-1:
                print ''

def OutputFormatTitle (Title, Nsl, Continue, num=0, map='multiple', format_='{0:4s} {1:5s}', word_add=' '):
        '''     Description:
        Function for adjusting formated output title with printed 
        output application data
        '''
        for number in range(Nsl):
            if number == 0:
                print Title
            print format_.format (str(number), word_add), 
        if map=='single':
            for number in range(int(num)):
                print format_.format (str(number), Continue), 
        elif map=='multiple':
            for txt in Continue:
                print '  %s  ' % txt,
        print '\n'
