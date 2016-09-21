import numpy as np
import matplotlib.pyplot as plt

N = 21
min = 10
max = 0
count = 0
count2 = 0
main_chain = np.zeros((N,), dtype=np.int)
side_chain = np.zeros((N,), dtype=np.int)
ring_chain = np.zeros((N,), dtype=np.int)
main_chain4 = np.zeros((N,), dtype=np.int)
side_chain4 = np.zeros((N,), dtype=np.int)
ring_chain4 = np.zeros((N,), dtype=np.int)
dict_residues = {'ALA':0, 'ARG':1, 'ASN':2, 'ASP':3, 'CYS':4, 'GLN':5, 'GLU':6, 'GLY':7, 'HIS':8, 'ILE':9, 'LEU':10, 'LYS':11, 'MET':12, 'PHE':13, 'PRO':14, 'SER':15, 'THR':16, 'TRP':17, 'TYR':18, 'VAL':19}
dict_main = {}
dict_side = {}
dict_ring = {}
n1 = 0
n2 = 0
n3 = 0
ring_id = 'a'
nA = 0
with open('new_sulfur.txt', 'r') as f:
    for line in f:
        newLine = ' '.join(line.split()).split(' ')
        repeat_main = False
        repeat_side = False
        repeat_ring = False
        if newLine[3] == newLine[7] and newLine[3] == newLine[11] and (not newLine[16] == 'DT') and (not newLine[16] == 'DA') and (not newLine[16] == 'DC') and (not newLine[16] == 'DG') and (not newLine[16] == 'A') and (not newLine[16] == 'G') and (not newLine[16] == 'C') and (not newLine[16] == 'U') and (not newLine[0] == 'PDB') and (not newLine[18] == 'H'):
            if newLine[15] == 'O' or newLine[15] == 'C' or newLine[15] == 'CA' or newLine[15] == 'N':
                for i in range(0,n1):
                    if dict_main[i]['id'] == newLine[0] and dict_main[i]['ligand_name'] == newLine[11] and dict_main[i]['res_id'] == newLine[17]:
                        repeat_main = True
                if not repeat_main:
                    dict_main[n1] = {'id':newLine[0],'ligand_name':newLine[11], 'res_id':newLine[17]}
                    n1+=1
                    if float(newLine[22]) <= 3.5:
                        main_chain[dict_residues[newLine[16]]+1]+=1
                    if float(newLine[22]) <= 4.0:
                        main_chain4[dict_residues[newLine[16]]+1]+=1
            elif newLine[14] == 'R':
                if not newLine[0] == ring_id:
                    nA+=1
                    ring_id = newLine[0]
                for i in range(0,n3):
                    if dict_ring[i]['id'] == newLine[0] and dict_ring[i]['ligand_name'] == newLine[11] and dict_ring[i]['res_id'] == newLine[16]:
                        repeat_ring = True
                if not repeat_ring:
                    dict_ring[n3] = {'id':newLine[0],'ligand_name':newLine[11], 'res_id':newLine[16]}
                    n3+=1
                    if float(newLine[19]) <= 5.0:
                        ring_chain[dict_residues[newLine[15]]+1]+=1
                    if float(newLine[19]) <= 6.0:
                        ring_chain4[dict_residues[newLine[15]]+1]+=1
            elif not newLine[14] == 'R':
                for i in range(0,n2):
                    if dict_side[i]['id'] == newLine[0] and dict_side[i]['ligand_name'] == newLine[11] and dict_side[i]['res_id'] == newLine[17]:
                        repeat_side = True
                if not repeat_side:
#                    if newLine[16] == 'PHE':
#                        print newLine[0]
#                    if newLine[16] == 'TYR':
#                        print '\t\t',newLine[0]
                    dict_side[n2] = {'id':newLine[0],'ligand_name':newLine[11], 'res_id':newLine[17]}
                    n2+=1
                    if float(newLine[22]) <= 3.5:
                        side_chain[dict_residues[newLine[16]]+1]+=1
                    if float(newLine[22]) <= 4.0:
                        side_chain4[dict_residues[newLine[16]]+1]+=1

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
labelx = -0.3

fig1 = plt.figure()
##fig, axarr = plt.subplots(3,3)
#
dist = main_chain
dist1 = main_chain4
ax0 = plt.subplot(121)
rects = ax0.bar(ind, dist, width, color='b')
rects1 = ax0.bar(ind+width, dist1, width, color='r')
ax0.legend( (rects[0], rects1[0]), ('Interactions $\mathregular{\leq}$ 3.5A', 'Interactions $\mathregular{\leq}$ 4.0A'), loc='upper left' )
ax0.set_xlim([1, N])
#
ax0.set_title('Main Chain Interactions', fontweight='bold', fontsize='20')
ax0.set_ylabel('No. of Interactions',fontweight='bold', fontsize='16')
ax0.set_xlabel('Amino Acids',fontweight='bold', fontsize='16')
for tick in ax0.yaxis.get_ticklabels():
    tick.set_weight('bold')
    tick.set_fontsize('14')
ax0.set_xticks(ind+width)
ax0.set_xticklabels(('', 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'), fontweight='bold', fontsize='14', rotation='vertical')
#plt.setp(ax0.get_xticklabels(), visible=False)
#plt.setp(ax0.get_xticklabels()[::-2], visible=True)
#xticks = ax0.xaxis.get_major_ticks()
#xticks[1].label1.set_visible(True)
#




#fig2 = plt.figure()
##fig, axarr = plt.subplots(3,3)
#
dist = side_chain
dist1 = side_chain4
ax1 = plt.subplot(122)
rects = ax1.bar(ind, dist, width, color='b')
rects1 = ax1.bar(ind+width, dist1, width, color='r')
ax1.legend( (rects[0], rects1[0]), ('Interactions $\mathregular{\leq}$ 3.5A', 'Interactions $\mathregular{\leq}$ 4.0A'), loc='best' )
ax1.set_xlim([1, N])
#
ax1.set_title('Side Chain Interactions', fontweight='bold', fontsize='20')
ax1.set_ylabel('No. of Interactions',fontweight='bold', fontsize='16')
ax1.set_xlabel('Amino Acids',fontweight='bold', fontsize='16')
for tick in ax1.yaxis.get_ticklabels():
    tick.set_weight('bold')
    tick.set_fontsize('14')
ax1.set_xticks(ind+width)
ax1.set_xticklabels(('','ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'), fontweight='bold', fontsize='14', rotation='vertical')
#plt.setp(ax1.get_xticklabels(), visible=False)
#plt.setp(ax1.get_xticklabels()[::-2], visible=True)

fig2 = plt.figure()
##fig, axarr = plt.subplots(3,3)
#
dist = ring_chain
dist1 = ring_chain4
ax0 = plt.subplot(111)
rects = ax0.bar(ind, dist, width, color='b')
rects1 = ax0.bar(ind+width, dist1, width, color='r')
ax0.legend( (rects[0], rects1[0]), ('Interactions $\mathregular{\leq}$ 5.0A', 'Interactions $\mathregular{\leq}$ 6.0A'), loc='upper left' )
ax0.set_xlim([1, N])
#
ax0.set_title('Ring Interactions', fontweight='bold', fontsize='20')
ax0.set_ylabel('No. of Interactions',fontweight='bold', fontsize='16')
ax0.set_xlabel('Amino Acids',fontweight='bold', fontsize='16')
for tick in ax0.yaxis.get_ticklabels():
    tick.set_weight('bold')
    tick.set_fontsize('14')
ax0.set_xticks(ind+width)
ax0.set_xticklabels(('', 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'), fontweight='bold', fontsize='14', rotation='vertical')


#fig3 = plt.figure()
###fig, axarr = plt.subplots(3,3)
##
#dist = Aromatic
#ax2 = plt.subplot(211)
#rects = ax2.bar(ind+width, dist, width, color='b')
#ax2.set_xlim([1, N])
##
#ax2.set_title('Aromatic', fontweight='bold', fontsize='20')
#ax2.set_ylabel('Number of Interactions',fontweight='bold', fontsize='16')
#ax2.set_xlabel('Bond Angle',fontweight='bold', fontsize='16')
#for tick in ax2.yaxis.get_ticklabels():
#    tick.set_weight('bold')
#    tick.set_fontsize('14')
#ax2.set_xticks(ind+width)
#ax2.set_xticklabels(('0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150','160','170','180'), fontweight='bold', fontsize='14')
##plt.setp(ax2.get_xticklabels(), visible=False)
##plt.setp(ax2.get_xticklabels()[::-5], visible=True)
#
#dist = Aromatic_pi
#ax3 = plt.subplot(212)
#rects = ax3.bar(ind+width, dist, width, color='b')
#ax3.set_xlim([1, N])
##
#ax3.set_title('Aromatic Pi', fontweight='bold', fontsize='20')
#ax3.set_ylabel('Number of Interactions',fontweight='bold', fontsize='16')
#ax3.set_xlabel('Bond Angle',fontweight='bold', fontsize='16')
#for tick in ax3.yaxis.get_ticklabels():
#    tick.set_weight('bold')
#    tick.set_fontsize('14')
#ax3.set_xticks(ind+width)
#ax3.set_xticklabels(('0','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90'), fontweight='bold', fontsize='14')
#plt.setp(ax2.get_xticklabels(), visible=False)
#plt.setp(ax2.get_xticklabels()[::-5], visible=True)


#
fig1.set_tight_layout(True)
print 'Number of side chain nteractions:',sum(side_chain)
print 'Number of main chain nteractions:',sum(main_chain)
print 'Number of ring nteractions:',sum(ring_chain4)
print 'Number of unique ring nteractions:',nA
#fig2.set_tight_layout(True)
#fig3.set_tight_layout(True)
plt.show()







