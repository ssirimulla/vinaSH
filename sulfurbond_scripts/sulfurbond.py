
#! /usr/bin/env python

import re
import math
import numpy as np
import itertools
import sys
import os
from Tkinter import Tk
from Tkinter import *
import ttk
from tkFileDialog import askdirectory

Tk().withdraw()


userpdb=askdirectory(title='Select PDB containing folder')
pdbfiles=os.listdir(userpdb)
useroutpath=askdirectory(title='Where would you like the ouput to be saved?')

useroutname=raw_input('Enter desired name for output file: ')+'.txt'

#the desired name and filepath of output file
rawoutfile=open(useroutpath+'/'+useroutname,'w')
rawoutfile.close()
rawoutfile=open(useroutpath+'/'+useroutname,'a')
cycnum=0


#outward loop
#repeats searh for each file in pdb folder
print("Checking pdb files...")
header='{:4s} {:6s} {:3s} {:3s} {:2s} {:6s} {:3s} {:3s} {:2s} {:6s} {:4s} {:3s} {:3s} {:1s} {:6s} {:3s} {:3s} {:4s} {:2s} {:6s} {:6s} {:6s} {:6s} {:6s}'.format('PDB ','Het#   ',
                                                                                                                            'Het','Res',
                                                                                                                            'Ty','Het#   ','Het','Res',
                                                                                                                            'Ty','Het#  ',
                                                                                                                            'Het ','Res',
                                                                                                                            'Ty ','|','Pro#   ',
                                                                                                                            'Atm','Res',
                                                                                                                            'R# ','Ty',
                                                                                                                            'Theta1','Theta2','Theta3','  Dis ','C-X--pi')
rawoutfile.write(header+ '\n')


for fil in pdbfiles:
    
    InFileName = userpdb +'/'+ fil
    
    InFile = open(InFileName, 'r')
    print InFileName[-8:-4]
    dict_hetsulfur = { }
    dict_hetsul={}
    dict_hetcarbon1= { }
    dict_hetcarbon2= { }
    dict_hetcarbon_holder= { }
    dict_conect = { }
    dict_atom = { }
    modres=[]
    
    i = 0
    for line in InFile:
            record_type=line[0:6]
            atom_type=line[76:78]
            newLine = ' '.join(line.split()).split(' ')
            if record_type == 'MODRES':
                modres.append(newLine[2])
            ## Get Sulfur information
            if record_type == 'HETATM' and (atom_type == ' S'):
                dict_hetsulfur[i] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20], 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78],'carbon':'False'}
                    #print dict_hetsulfur[i]
                    #print line
                i = i + 1
    f = i
    #print (dict_hetsulfur)
    InFile.close()

    InFile = open(InFileName, 'r')
    
    #Benzene Dictionaries

    dict_tyr = { }
    dict_phe = { }
    dict_trp = { }
    dict_cen_benz={ }
    
    ty=0
    ph=0
    tr=0
    ## Get ring information
    for line in InFile:
            record_type=line[0:6]
            res_name=line[17:20]
            side_chain=line[12:16]
            if record_type == 'ATOM  ' and res_name == 'TYR' and (side_chain==' CG 'or side_chain==' CD1' or side_chain==' CD2' or side_chain==' CE1' or side_chain==' CE2' or side_chain==' CZ '):
                dict_tyr[ty] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20],
                    'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                ty += 1
            elif record_type == 'ATOM  ' and res_name == 'PHE' and (side_chain==' CG 'or side_chain==' CD1' or side_chain==' CD2' or side_chain==' CE1' or side_chain==' CE2' or side_chain==' CZ '):
                dict_phe[ph] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20],
                    'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                ph += 1
            elif record_type == 'ATOM  ' and res_name == 'TRP' and (side_chain==' CD2'or side_chain==' CE2' or side_chain==' CE3' or side_chain==' CZ2' or side_chain==' CZ3' or side_chain==' CH2'):
                dict_trp[tr] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20],
                    'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                tr += 1

    InFile.close()
    
    #Finding the center and normal vector of benzene rings
    def benz_centerx(x1,x2,x3,x4,x5,x6):
        return (x1+x2+x3+x4+x5+x6)/6
    def benz_centery(y1,y2,y3,y4,y5,y6):
        return (y1+y2+y3+y4+y5+y6)/6
    def benz_centerz(z1,z2,z3,z4,z5,z6):
        return (z1+z2+z3+z4+z5+z6)/6
    
    
    def normal_v(x1,x2,x3,y1,y2,y3,z1,z2,z3):
        v1=[x2-x1,y2-y1,z2-z1]
        v2=[x3-x1,y3-y1,z3-z1]
        return np.cross(v1,v2)

#dictionary with only one index. used for loops
    dict_cycle={}
    
    benzall=0

    if ty%6 ==0:
        ty+=1
        for carb in range(0,ty//6):
            dict_cycle[0]={'cenx': benz_centerx(dict_tyr[6*carb]['x'],dict_tyr[(6*carb)+1]['x'],dict_tyr[(6*carb)+2]['x'],dict_tyr[(6*carb)+3]['x'],dict_tyr[(6*carb)+4]['x'],dict_tyr[(6*carb)+5]['x']),
            'ceny': benz_centery(dict_tyr[6*carb]['y'],dict_tyr[(6*carb)+1]['y'],dict_tyr[(6*carb)+2]['y'],dict_tyr[(6*carb)+3]['y'],dict_tyr[(6*carb)+4]['y'],dict_tyr[(6*carb)+5]['y']),
                'cenz': benz_centerz(dict_tyr[6*carb]['z'],dict_tyr[(6*carb)+1]['z'],dict_tyr[(6*carb)+2]['z'],dict_tyr[(6*carb)+3]['z'],dict_tyr[(6*carb)+4]['z'],dict_tyr[(6*carb)+5]['z'])
                }
            dict_cen_benz[benzall]={'cenx': dict_cycle[0]['cenx'],
                'ceny': dict_cycle[0]['ceny'] ,
                    'cenz': dict_cycle[0]['cenz'],
                        'norm_v': normal_v(dict_cycle[0]['cenx'],dict_tyr[6*carb]['x'],dict_tyr[(6*carb)+1]['x'],dict_cycle[0]['ceny'],dict_tyr[carb*6]['y'],dict_tyr[(carb*6)+1]['y'],
                                           dict_cycle[0]['cenz'],dict_tyr[carb*6]['z'],dict_tyr[(carb*6)+1]['z']),
                            'resseq':dict_tyr[6*carb]['resseq'],
                                'resname':dict_tyr[6*carb]['resname'],
                                    'sn_range': [dict_tyr[6*carb]['serial_number'],dict_tyr[(6*carb)+5]['serial_number']]}
            
            
            benzall+=1
    if ph%6 ==0:
        ph+=1
        for carb in range(0,ph//6):
            dict_cycle[0]={'cenx': benz_centerx(dict_phe[6*carb]['x'],dict_phe[(6*carb)+1]['x'],dict_phe[(6*carb)+2]['x'],dict_phe[(6*carb)+3]['x'],dict_phe[(6*carb)+4]['x'],dict_phe[(6*carb)+5]['x']),
        'ceny': benz_centery(dict_phe[6*carb]['y'],dict_phe[(6*carb)+1]['y'],dict_phe[(6*carb)+2]['y'],dict_phe[(6*carb)+3]['y'],dict_phe[(6*carb)+4]['y'],dict_phe[(6*carb)+5]['y']),
            'cenz': benz_centerz(dict_phe[6*carb]['z'],dict_phe[(6*carb)+1]['z'],dict_phe[(6*carb)+2]['z'],dict_phe[(6*carb)+3]['z'],dict_phe[(6*carb)+4]['z'],dict_phe[(6*carb)+5]['z'])
            }
            dict_cen_benz[benzall]={'cenx': dict_cycle[0]['cenx'],
                'ceny': dict_cycle[0]['ceny'] ,
                    'cenz': dict_cycle[0]['cenz'],
                        'norm_v': normal_v(dict_cycle[0]['cenx'],dict_phe[6*carb]['x'],dict_phe[(6*carb)+1]['x'],dict_cycle[0]['ceny'],dict_phe[carb*6]['y'],dict_phe[(carb*6)+1]['y'],
                                           dict_cycle[0]['cenz'],dict_phe[carb*6]['z'],dict_phe[(carb*6)+1]['z']),
                            'resseq':dict_phe[6*carb]['resseq'],
                                'resname':dict_phe[6*carb]['resname'],
                                    'sn_range': [dict_phe[6*carb]['serial_number'],dict_phe[(6*carb)+5]['serial_number']]}
    
    
            benzall+=1
    if tr%6 ==0:
        tr+=1
        for carb in range(0,tr//6):
            dict_cycle[0]={'cenx': benz_centerx(dict_trp[6*carb]['x'],dict_trp[(6*carb)+1]['x'],dict_trp[(6*carb)+2]['x'],dict_trp[(6*carb)+3]['x'],dict_trp[(6*carb)+4]['x'],dict_trp[(6*carb)+5]['x']),
            'ceny': benz_centery(dict_trp[6*carb]['y'],dict_trp[(6*carb)+1]['y'],dict_trp[(6*carb)+2]['y'],dict_trp[(6*carb)+3]['y'],dict_trp[(6*carb)+4]['y'],dict_trp[(6*carb)+5]['y']),
                'cenz': benz_centerz(dict_trp[6*carb]['z'],dict_trp[(6*carb)+1]['z'],dict_trp[(6*carb)+2]['z'],dict_trp[(6*carb)+3]['z'],dict_trp[(6*carb)+4]['z'],dict_trp[(6*carb)+5]['z'])
                }
            dict_cen_benz[benzall]={'cenx': dict_cycle[0]['cenx'],
                'ceny': dict_cycle[0]['ceny'] ,
                    'cenz': dict_cycle[0]['cenz'],
                        'norm_v': normal_v(dict_cycle[0]['cenx'],dict_trp[6*carb]['x'],dict_trp[(6*carb)+1]['x'],dict_cycle[0]['ceny'],dict_trp[carb*6]['y'],dict_trp[(carb*6)+1]['y'],
                                           dict_cycle[0]['cenz'],dict_trp[carb*6]['z'],dict_trp[(carb*6)+1]['z']),
                            'resseq':dict_trp[6*carb]['resseq'],
                                'resname':dict_trp[6*carb]['resname'],
                                    'sn_range': [dict_trp[6*carb]['serial_number'],dict_trp[(6*carb)+5]['serial_number']]}
            
            
            benzall+=1


    InFile = open(InFileName, 'r')
    n=0
    ## Get connections to sulfur
    for line in InFile:
        record_type=line[0:6]
        serial_number=line[6:11]
        if record_type == 'CONECT':
            for q in dict_hetsulfur:
                if serial_number==dict_hetsulfur[q]['serial_number']:
                    dict_conect[n] = {'serial_number1':line[6:11], 'serial_number2':line[11:16],'serial_number3':line[16:21],'serial_number4':line[21:26]}
                    n = n + 1

    InFile.close()
    i = 0
    n= 0
    m=0
    InFile = open(InFileName, 'r')
    for line in InFile:
        record_type=line[0:6]
        serial_number=line[6:11]
        atom_type=line[76:78]
        res_name=line[17:20]
        found=False
        if record_type == 'HETATM' and atom_type == ' C' and (res_name not in modres):
            for con in dict_conect:
                ## Find a Sulfur that is connected to exactly two carbons
                if (serial_number==dict_conect[con]['serial_number2'] or serial_number==dict_conect[con]['serial_number3']) and dict_conect[con]['serial_number4']=='     ':
                    for sul in range(len(dict_hetsulfur)):
                        ## Second Carbon
                        if dict_hetsulfur[sul]['serial_number']==dict_conect[con]['serial_number1'] and dict_hetsulfur[sul]['carbon']=='One' and serial_number==dict_conect[con]['serial_number3']:
                            dict_hetsul[n]=dict_hetsulfur[sul]
                            print dict_conect[con]
                            dict_hetsulfur[sul]['carbon'] = 'True'
                            print dict_hetsulfur[sul]['serial_number']
                            
                            for j in range(0,m):
                                if dict_hetcarbon_holder[j]['sul_serial'] == dict_hetsulfur[sul]['serial_number'] and serial_number != dict_hetcarbon_holder[j]['serial_number']:
                                    dict_hetcarbon1[n] = dict_hetcarbon_holder[j]
                                    dict_hetcarbon2[n] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17],
                                        'resname':line[17:20], 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27],
                                            'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78],
                                                'H1':None,'H2':None,'H3':None}
                                    n = n + 1
                                    break
                            break
                        
                        ## First Carbon
                        elif dict_hetsulfur[sul]['serial_number']==dict_conect[con]['serial_number1'] and dict_hetsulfur[sul]['carbon']=='False' and serial_number==dict_conect[con]['serial_number2']:
                            dict_hetsulfur[sul]['carbon'] = 'One'
    
                            dict_hetcarbon_holder[m] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17],
                                'resname':line[17:20], 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27],
                                    'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78],
                                        'H1':None,'H2':None,'H3':None,'sul_serial':dict_hetsulfur[sul]['serial_number']}
                            m = m + 1
                            break
    w = n
    InFile.close()

    i = 0
    d = 0
    s = 0
    ## Angles with non-rings
    for j in range(0,w):
            InFile = open(InFileName, 'r')
            for line in InFile:
                    record_type=line[0:6]
                    if record_type == 'ATOM  ' and (dict_hetsul[j]['resname'] not in modres):
                        try:
                            distance = math.sqrt(((dict_hetsul[j]['x'] - float(line[30:38]))**2) + ((dict_hetsul[j]['y'] - float(line[38:46]))**2) + ((dict_hetsul[j]['z'] - float(line[46:54]))**2)) # calculate length between sulfur and atoms in protein
                        except(ZeroDivisionError,ValueError):
                                continue
                        if distance < 4.0: # if length is less than 4.0
                            dict_atom[s] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20], 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}

                            v1 = [(dict_hetsul[j]['x']-dict_hetcarbon1[j]['x']),(dict_hetsul[j]['y']-dict_hetcarbon1[j]['y']),(dict_hetsul[j]['z']-dict_hetcarbon1[j]['z'])] # vector v1 (hetsulfur - hetcarbon1)
                            v2 = [(dict_hetsul[j]['x']-dict_atom[s]['x']),(dict_hetsul[j]['y']-dict_atom[s]['y']),(dict_hetsul[j]['z']-dict_atom[s]['z'])] # vector v2 (hetsulfur- atom)
                            
                            v3 = [(dict_hetsul[j]['x']-dict_hetcarbon2[j]['x']),(dict_hetsul[j]['y']-dict_hetcarbon2[j]['y']),(dict_hetsul[j]['z']-dict_hetcarbon2[j]['z'])] # vector v1 (hetsulfur - hetcarbon2)
                            def dotproduct(v1, v2):
                                        return sum((a*b) for a, b in zip(v1, v2))
                            def dotproduct2(v3, v2):
                                return sum((a*b) for a, b in zip(v3, v2))
                            def dotproduct3(v1, v3):
                                return sum((a*b) for a, b in zip(v1, v3))
                            def length(v):
                                        return math.sqrt(dotproduct(v, v))
                                    
                            length_v1 = length(v1) #calculate length of vector 1
                            length_v2 = length(v2) #calculate length of vector 2
                            length_v3 = length(v3) #calculate length of vector 3
                                    
                            def angle1(v1, v2):
                                        return math.acos(dotproduct(v1, v2) / (length_v1 * length_v2))
                            def angle2(v3, v2):
                                return math.acos(dotproduct2(v3, v2) / (length_v3 * length_v2))
                            def angle3(v1, v3):
                                return math.acos(dotproduct3(v1, v3) / (length_v1 * length_v3))
                            try:
                                        Angle1 = angle1(v1, v2) *180/3.141592
                            except(ZeroDivisionError,ValueError):
                                        continue
                            try:
                                Angle2 = angle2(v3, v2) *180/3.141592
                            except(ZeroDivisionError,ValueError):
                                continue
                            try:
                                Angle3 = angle3(v1, v3) *180/3.141592
                            except(ZeroDivisionError,ValueError):
                                continue
                            if (Angle3 <= 180 and Angle3 >= 0):
                                rawout='{:4s} {:6d} {:3s} {:3s} {:2s} {:6d} {:3s} {:3s} {:2s} {:6d} {:4s} {:3s} {:3s} {:1s} {:6d} {:3s} {:3s} {:4d} {:2s} {:6.2f} {:6.2f} {:6.2f} {:6.2f}'.format(InFileName[-8:-4],int(dict_hetcarbon1[j]['serial_number']),
                                                                                                                                                              dict_hetcarbon1[j]['fullname'],dict_hetcarbon1[j]['resname'],
                                                                                                                                                                                dict_hetcarbon1[j]['atom_type'],int (dict_hetcarbon2[j]['serial_number']), dict_hetcarbon2[j]['fullname'],dict_hetcarbon2[j]['resname'],
                                                                                                                                                                                dict_hetcarbon2[j]['atom_type'], int(dict_hetsul[j]['serial_number']),
                                                                                                                                                              dict_hetsul[j]['fullname'],dict_hetsul[j]['resname'],
                                                                                                                                                              dict_hetsul[j]['atom_type'],'|',int(dict_atom[s]['serial_number']),
                                                                                                                                                              dict_atom[s]['fullname'],dict_atom[s]['resname'],
                                                                                                                                                              int(dict_atom[s]['resseq']),dict_atom[s]['atom_type'],
                                                                                                                                                              round(Angle1,2),round(Angle2,2),round(Angle3,2),round(distance,2))
                                rawoutfile.write(rawout+'\n')


    InFile.close()
    
    #Finding benzine rings within 6a
    rings=True
    if benzall==0:
        rings=False

    def benz_dis(x1,x2,y1,y2,z1,z2):
        return math.sqrt((x2-x1)**2 +(y2-y1)**2 +(z2-z1)**2)
    
    def plane_angle(x1,x2,y1,y2,z1,z2,normp):
        hv=[x2-x1,y2-y1,z2-z1]
        return round(math.degrees(math.asin((np.dot(hv,normp))/(benz_dis(hv[0],0,hv[1],0,hv[2],0)*benz_dis(normp[0],0,normp[1],0,normp[2],0)))),2)
    def bh_bond_angle1(hv,bv):
        return round(math.degrees(math.acos((np.dot(hv, bv) / ((benz_dis(hv[0],0,hv[1],0,hv[2],0) * benz_dis(bv[0],0,bv[1],0,bv[2],0)))))),2)
    def bh_bond_angle2(hv1,bv1):
        return round(math.degrees(math.acos((np.dot(hv1, bv1) / ((benz_dis(hv1[0],0,hv1[1],0,hv1[2],0) * benz_dis(bv1[0],0,bv1[1],0,bv1[2],0)))))),2)
    def csc_bond_angle(v1,v2):
        return round(math.degrees(math.acos((np.dot(v1, v2) / ((benz_dis(v1[0],0,v1[1],0,v1[2],0) * benz_dis(v2[0],0,v2[1],0,v2[2],0)))))),2)
    
    ## Angles with rings
    for ent in range(0,w):
        for b in range(0,benzall):
            try:
                bdis=benz_dis(dict_cen_benz[b]['cenx'],dict_hetsul[ent]['x'],
                              dict_cen_benz[b]['ceny'],dict_hetsul[ent]['y'],
                              dict_cen_benz[b]['cenz'],dict_hetsul[ent]['z'])
                    
                bpangle=plane_angle(dict_cen_benz[b]['cenx'],dict_hetsul[ent]['x'],
                                    dict_cen_benz[b]['ceny'],dict_hetsul[ent]['y'],
                                    dict_cen_benz[b]['cenz'],dict_hetsul[ent]['z'],
                                    dict_cen_benz[b]['norm_v'])
                baangle1=bh_bond_angle1([dict_hetcarbon1[ent]['x']-dict_hetsul[ent]['x'],
                                         dict_hetcarbon1[ent]['y']-dict_hetsul[ent]['y'],
                                         dict_hetcarbon1[ent]['z']-dict_hetsul[ent]['z']],
                                        [dict_cen_benz[b]['cenx']-dict_hetsul[ent]['x'],
                                         dict_cen_benz[b]['ceny']-dict_hetsul[ent]['y'],
                                         dict_cen_benz[b]['cenz']-dict_hetsul[ent]['z']])
                baangle2=bh_bond_angle2([dict_hetcarbon2[ent]['x']-dict_hetsul[ent]['x'],
                                         dict_hetcarbon2[ent]['y']-dict_hetsul[ent]['y'],
                                         dict_hetcarbon2[ent]['z']-dict_hetsul[ent]['z']],
                                        [dict_cen_benz[b]['cenx']-dict_hetsul[ent]['x'],
                                         dict_cen_benz[b]['ceny']-dict_hetsul[ent]['y'],
                                         dict_cen_benz[b]['cenz']-dict_hetsul[ent]['z']])
                cscngle=csc_bond_angle([dict_hetcarbon1[ent]['x']-dict_hetsul[ent]['x'],
                                         dict_hetcarbon1[ent]['y']-dict_hetsul[ent]['y'],
                                         dict_hetcarbon1[ent]['z']-dict_hetsul[ent]['z']],
                                        [dict_hetcarbon2[ent]['x']-dict_hetsul[ent]['x'],
                                         dict_hetcarbon2[ent]['y']-dict_hetsul[ent]['y'],
                                         dict_hetcarbon2[ent]['z']-dict_hetsul[ent]['z']])
            except(ZeroDivisionError,ValueError):
                continue
            if rings and bdis<= 6: #and cscngle >= 105 and cscngle <= 110:
                rawout='{:4s} {:6d} {:3s} {:3s} {:2s} {:6d} {:3s} {:3s} {:2s} {:6d} {:4s} {:3s} {:3s} {:1s} {:6s} {:3s} {:3s} {:4d} {:2s} {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f}'.format(InFileName[-8:-4],int(dict_hetcarbon1[ent]['serial_number']),
                                                                                                                                                  dict_hetcarbon1[ent]['fullname'],dict_hetcarbon1[ent]['resname'],
                                                                                                                                                                            dict_hetcarbon1[ent]['atom_type'],int(dict_hetcarbon2[ent]['serial_number']),dict_hetcarbon2[ent]['fullname'],dict_hetcarbon2[ent]['resname'],
                                                                                                                                                                            dict_hetcarbon2[ent]['atom_type'],int(dict_hetsul[ent]['serial_number']),
                                                                                                                                                  dict_hetsul[ent]['fullname'],dict_hetsul[ent]['resname'],
                                                                                                                                                  dict_hetsul[ent]['atom_type'],'|','      ',
                                                                                                                                                  '  R ',dict_cen_benz[b]['resname'],
                                                                                                                                                  int(dict_cen_benz[b]['resseq']),' ',
                                                                                                                                                  round(baangle1,2),round(baangle2,2),round(bdis,2),round(bpangle,2),round(cscngle,2))
                rawoutfile.write(rawout+'\n')

rawoutfile.close() 
print('Sulfur interaction Analysis Complete')
