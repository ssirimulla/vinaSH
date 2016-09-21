import numpy as np
import matplotlib.pyplot as plt

N = 19
N1 = 181
min = 10
max = 0
count = 0
count2 = 0
Oxygen = np.zeros((N,), dtype=np.int)
Nitrogen = np.zeros((N,), dtype=np.int)
Oxygen_DNA = np.zeros((N,), dtype=np.int)
Nitrogen_DNA = np.zeros((N,), dtype=np.int)
Aromatic = np.zeros((N,), dtype=np.int)
Aromatic_pi = np.zeros((N,), dtype=np.int)
Oxygen1 = np.zeros((N1,), dtype=np.int)
Nitrogen1 = np.zeros((N1,), dtype=np.int)
#Oxygen_DNA1 = np.zeros((N1,), dtype=np.int)
#Nitrogen_DNA1 = np.zeros((N1,), dtype=np.int)
dict_O = {}
dict_N = {}
dict_A = {}
sulfur_id = 's'
oxygen_id = 'o'
nitrogen_id = 'n'
n1 = 0
n2 = 0
n3 = 0
nS = 0
nO = 0
nN = 0
with open('sulfur_no_angle_limit.txt', 'r') as f:
    for line in f:
        newLine = ' '.join(line.split()).split(' ')
        repeat_O = False
        repeat_N = False
        repeat_A = False
        if not newLine[0] == sulfur_id and newLine[3] == newLine[7] and newLine[3] == newLine[11]:
            nS+=1
            sulfur_id = newLine[0]
        if (not newLine[21] == 'Theta3') and float(newLine[21]) >= 0 and float(newLine[21]) < 181:
            if newLine[3] == newLine[7] and newLine[3] == newLine[11] and (not newLine[16] == 'DT') and (not newLine[16] == 'DA') and (not newLine[16] == 'DC') and (not newLine[16] == 'DG') and (not newLine[14] == 'R'):
                if newLine[18] == 'O' and float(newLine[22]) <= 3.32:
                    if not newLine[0] == oxygen_id:
                        nO+=1
                        oxygen_id = newLine[0]
                    for i in range(0,n1):
                        if dict_O[i]['id'] == newLine[0] and dict_O[i]['serial_number1'] == newLine[1] and dict_O[i]['serial_number2'] == newLine[5] and dict_O[i]['serial_number3'] == newLine[9] and dict_O[i]['serial_number4'] == newLine[14]:
                            repeat_O = True
                    if not repeat_O:
#                        print line
                        dict_O[n1] = {'id':newLine[0],'serial_number1':newLine[1], 'serial_number2':newLine[5], 'serial_number3':newLine[9], 'serial_number4':newLine[14]}
                        n1+=1
                        Oxygen1[int(round(float(newLine[19])))] += 1
                        Oxygen1[int(round(float(newLine[20])))] += 1
                        if float(newLine[19]) > 180 or float(newLine[20]) > 180:
                            print line
    #                    if float(newLine[19]) > 170 or float(newLine[20]) > 170:
    #                        print line
                        for i in range(0,N):
                            if float(newLine[19]) > i*10 and float(newLine[19]) <= i*10+10: Oxygen[i]+=1
                            if float(newLine[20]) > i*10 and float(newLine[20]) <= i*10+10: Oxygen[i]+=1
                elif newLine[18] == 'N' and float(newLine[22]) <= 3.35:
                    if newLine[18] == 'N' and not newLine[0] == nitrogen_id:
                        nN+=1
                        nitrogen_id = newLine[0]
                    for i in range(0,n2):
                        if dict_N[i]['id'] == newLine[0] and dict_N[i]['serial_number1'] == newLine[1] and dict_N[i]['serial_number2'] == newLine[5] and dict_N[i]['serial_number3'] == newLine[9] and dict_N[i]['serial_number4'] == newLine[14]:
                            repeat_N = True
                    if not repeat_N:
#                        print line
                        dict_N[n2] = {'id':newLine[0],'serial_number1':newLine[1], 'serial_number2':newLine[5], 'serial_number3':newLine[9], 'serial_number4':newLine[14]}
                        n2+=1
                        Nitrogen1[int(round(float(newLine[19])))] += 1
                        Nitrogen1[int(round(float(newLine[20])))] += 1
    #                    if float(newLine[19]) > 170 or float(newLine[20]) > 170:
    #                        print line
                        for i in range(0,N):
                            if float(newLine[19]) > i*10 and float(newLine[19]) <= i*10+10: Nitrogen[i]+=1
                            if float(newLine[20]) > i*10 and float(newLine[20]) <= i*10+10: Nitrogen[i]+=1
            elif newLine[3] == newLine[7] and newLine[3] == newLine[11] and (not newLine[14] == 'R'):
                if newLine[18] == 'O' and float(newLine[22]) <= 3.32:
                    for i in range(0,n1):
                        if dict_O[i]['id'] == newLine[0] and dict_O[i]['serial_number1'] == newLine[1] and dict_O[i]['serial_number2'] == newLine[5] and dict_O[i]['serial_number3'] == newLine[9] and dict_O[i]['serial_number4'] == newLine[14]:
                            repeat_O = True
                    if not repeat_O:
#                        if newLine[0] == '108d' and (newLine[1] == '513' or newLine[1] == '542'):
#                        print line
                        dict_O[n1] = {'id':newLine[0],'serial_number1':newLine[1], 'serial_number2':newLine[5], 'serial_number3':newLine[9], 'serial_number4':newLine[14]}
                        n1+=1
                        for i in range(0,N):
                            if float(newLine[19]) > i*10 and float(newLine[19]) <= i*10+10: Oxygen_DNA[i]+=1
                            if float(newLine[20]) > i*10 and float(newLine[20]) <= i*10+10: Oxygen_DNA[i]+=1
                elif newLine[18] == 'N' and float(newLine[22]) <= 3.35:
                    for i in range(0,n2):
                        if dict_N[i]['id'] == newLine[0] and dict_N[i]['serial_number1'] == newLine[1] and dict_N[i]['serial_number2'] == newLine[5] and dict_N[i]['serial_number3'] == newLine[9] and dict_N[i]['serial_number4'] == newLine[14]:
                            repeat_N = True
                    if not repeat_N:
#                        print line
                        dict_N[n2] = {'id':newLine[0],'serial_number1':newLine[1], 'serial_number2':newLine[5], 'serial_number3':newLine[9], 'serial_number4':newLine[14]}
                        n2+=1
                        for i in range(0,N):
                            if float(newLine[19]) > i*10 and float(newLine[19]) <= i*10+10: Nitrogen_DNA[i]+=1
                            if float(newLine[20]) > i*10 and float(newLine[20]) <= i*10+10: Nitrogen_DNA[i]+=1
            elif newLine[3] == newLine[7] and newLine[3] == newLine[11] and newLine[14] == 'R':
                for i in range(0,n3):
                    if dict_A[i]['id'] == newLine[0] and dict_A[i]['serial_number1'] == newLine[1] and dict_A[i]['serial_number2'] == newLine[5] and dict_A[i]['serial_number3'] == newLine[9] and dict_A[i]['aromatic'] == newLine[15]:
                        repeat_A = True
                if not repeat_A:
                    dict_A[n3] = {'id':newLine[0],'serial_number1':newLine[1], 'serial_number2':newLine[5], 'serial_number3':newLine[9], 'aromatic':newLine[15]}
                    n3+=1
                    for i in range(0,N):
    #                    if float(newLine[17]) > 170 or float(newLine[18]) > 170:
    #                        print line
                        if float(newLine[17]) > i*10 and float(newLine[17]) <= i*10+10: Aromatic[i]+=1
                        if float(newLine[18]) > i*10 and float(newLine[18]) <= i*10+10: Aromatic[i]+=1
                        if abs(float(newLine[20])) > i*5 and abs(float(newLine[20])) <= i*5+5: Aromatic_pi[i]+=1


#for i in range(0,n1):
#    print dict_O[i]

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
labelx = -0.3

fig1 = plt.figure(figsize=(19,7))
##fig, axarr = plt.subplots(3,3)
#
dist = Oxygen
dist1 = Oxygen_DNA
ax0 = plt.subplot(121)
rects = ax0.bar(ind, dist, width, color='b')
rects1 = ax0.bar(ind+width, dist1, width, color='r')
ax0.legend( (rects[0], rects1[0]), ('Protein interactions', 'DNA interactions'), loc='upper left' )
ax0.set_xlim([1, N])
#
ax0.set_title('S$\cdots$O', fontweight='bold', fontsize='20')
ax0.set_ylabel('No. of Interactions',fontweight='bold', fontsize='16')
ax0.set_xlabel('Bond Angle',fontweight='bold', fontsize='16')
for tick in ax0.yaxis.get_ticklabels():
    tick.set_weight('bold')
    tick.set_fontsize('14')
ax0.set_xticks(ind+width)
ax0.set_xticklabels(('0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150','160','170','180'), fontweight='bold', fontsize='14')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax0.get_xticklabels()[::-2], visible=True)
#xticks = ax0.xaxis.get_major_ticks()
#xticks[1].label1.set_visible(True)
#
#
dist = Nitrogen
dist1 = Nitrogen_DNA
ax1 = plt.subplot(122)
rects = ax1.bar(ind, dist, width, color='b')
rects1 = ax1.bar(ind+width, dist1, width, color='r')
ax1.legend( (rects[0], rects1[0]), ('Protein interactions', 'DNA interactions'), loc='upper left' )
ax1.set_xlim([1, N])
#
ax1.set_title('S$\cdots$N', fontweight='bold', fontsize='20')
ax1.set_ylabel('No. of Interactions',fontweight='bold', fontsize='16')
ax1.set_xlabel('Bond Angle',fontweight='bold', fontsize='16')
for tick in ax1.yaxis.get_ticklabels():
    tick.set_weight('bold')
    tick.set_fontsize('14')
ax1.set_xticks(ind+width)
ax1.set_xticklabels(('0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150','160','170','180'), fontweight='bold', fontsize='14')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels()[::-2], visible=True)

#plt.savefig('Atomic_level_sulfur_interactions_with_angle_restriction_125_130.tiff', format='tiff', dpi=300)
tick_lables = np.linspace(0, 180, num=N1)
tick_lables = tick_lables.astype(int)
#print tick_lables
fig2 = plt.figure(figsize=(19,7))
##fig, axarr = plt.subplots(3,3)
#
ind = np.arange(N1)
dist = np.cumsum(Oxygen1/float(np.sum(Oxygen1)))
ax0 = plt.subplot(121)
rects = ax0.plot(ind, dist, linewidth=3, color='b')
#ax0.legend( (rects[0], rects1[0]), ('Protein interactions', 'DNA interactions'), loc='upper left' )
ax0.set_xlim([1, N1])
#
ax0.set_title('S$\cdots$O', fontweight='bold', fontsize='20')
#ax0.set_ylabel('No. of Interactions',fontweight='bold', fontsize='16')
ax0.set_xlabel('Bond Angle',fontweight='bold', fontsize='16')
for tick in ax0.yaxis.get_ticklabels():
    tick.set_weight('bold')
    tick.set_fontsize('14')
ax0.set_xticks(ind+width)
ax0.set_xticklabels(tick_lables, fontweight='bold', fontsize='14')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax0.get_xticklabels()[::-10], visible=True)
#xticks = ax0.xaxis.get_major_ticks()
#xticks[1].label1.set_visible(True)

#
dist = np.cumsum(Nitrogen1/float(np.sum(Nitrogen1)))
#dist1 = Nitrogen_DNA
ax1 = plt.subplot(122)
rects = ax1.plot(ind, dist, linewidth=3, color='b')
#rects1 = ax1.bar(ind+width, dist1, width, color='r')
#ax1.legend( (rects[0], rects1[0]), ('Protein interactions', 'DNA interactions'), loc='upper left' )
ax1.set_xlim([1, N1])
ax1.set_ylim([0, 1])
#
ax1.set_title('S$\cdots$N', fontweight='bold', fontsize='20')
#ax1.set_ylabel('No. of Interactions',fontweight='bold', fontsize='16')
ax1.set_xlabel('Bond Angle',fontweight='bold', fontsize='16')
for tick in ax1.yaxis.get_ticklabels():
    tick.set_weight('bold')
    tick.set_fontsize('14')
ax1.set_xticks(ind+width)
ax1.set_xticklabels(tick_lables, fontweight='bold', fontsize='14')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels()[::-10], visible=True)


fig3 = plt.figure(figsize=(17,7))
ax0 = plt.subplot(111)
#
ind = np.arange(N1)
dist = np.cumsum(Oxygen1/float(np.sum(Oxygen1)))
dist1 = np.cumsum(Nitrogen1/float(np.sum(Nitrogen1)))
rects = ax0.plot(ind, dist, linewidth=3, color='b')
rects1 = ax0.plot(ind, dist1, linewidth=3, color='r')
ax0.legend( (rects[0], rects1[0]), ('S$\cdots$O', 'S$\cdots$N'), loc='upper left' )
ax0.set_xlim([1, N1])
ax0.set_ylim([0, 1])

ax0.set_xlabel('Bond Angle',fontweight='bold', fontsize='16')
for tick in ax0.yaxis.get_ticklabels():
    tick.set_weight('bold')
    tick.set_fontsize('14')
ax0.set_xticks(ind+width)
ax0.set_xticklabels(tick_lables, fontweight='bold', fontsize='14')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax0.get_xticklabels()[::-10], visible=True)
#xticks = ax0.xaxis.get_major_ticks()
#xticks[1].label1.set_visible(True)

#


#plt.savefig('Old_5_Sulfur_JCIM_CDF.tiff', format='tiff', dpi=300)
#print n1, n2


#print 'Number of Oxygen bonds:', sum(Oxygen)+sum(Oxygen_DNA)
#print 'Number of Nitrogen bonds:', sum(Nitrogen)+sum(Nitrogen_DNA)
#print 'Number of unique Oxygen interactions:', nO
#print 'Number of unique Nitrogen interactions:', nN
#print 'Number of Sulfinated ligands:', nS

#
fig1.set_tight_layout(True)
fig2.set_tight_layout(True)
fig3.set_tight_layout(True)
plt.show()







