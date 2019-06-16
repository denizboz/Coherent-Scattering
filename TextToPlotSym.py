#This program does the same job with TextToPlot.py, only for symmetric combinations of nuclear spin arrangements.

import numpy as np
from matplotlib import pyplot as plt

Pref0011 = np.loadtxt("Nuc4-0011.txt")
Pref0101 = np.loadtxt("Nuc4-0101.txt")
Pref1010 = np.loadtxt("Nuc4-1010.txt")
Pref1100 = np.loadtxt("Nuc4-1100.txt")

Pref000111 = np.loadtxt("Nuc6-000111.txt")
Pref010101 = np.loadtxt("Nuc6-010101.txt")
Pref101010 = np.loadtxt("Nuc6-101010.txt")
Pref111000 = np.loadtxt("Nuc6-111000.txt")

Pref00001111 = np.loadtxt("Nuc8-00001111.txt")
Pref01010101 = np.loadtxt("Nuc8-01010101.txt")
Pref10101010 = np.loadtxt("Nuc8-10101010.txt")
Pref11110000 = np.loadtxt("Nuc8-11110000.txt")

r = np.linspace(0,1,101)
rsq = r**2

plt.figure(1)
plt.subplots_adjust(left=.1, bottom=.11, right=.97, top=.97, wspace=.3, hspace=.4)

plt.subplot(2,2,1)
plt.plot(rsq,Pref0011)
plt.plot(rsq,Pref000111)
plt.plot(rsq,Pref00001111)
plt.legend([r'$\downarrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\uparrow$',
r'$\downarrow\!\!\!\!\downarrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\uparrow\!\!\!\!\uparrow$',
r'$\downarrow\!\!\!\!\downarrow\!\!\!\!\downarrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\uparrow\!\!\!\!\uparrow\!\!\!\!\uparrow$'],fontsize=8)
#plt.title('00001111',fontsize=14)
plt.xlabel("$|r|^2$",fontsize=12)
plt.ylabel("$P_r$",fontsize=12)
plt.grid(b=None, which='major', axis='both')
plt.xlim((0,1))
plt.ylim((0,1.1))
#plt.tight_layout()

plt.subplot(2,2,2)
plt.plot(rsq,Pref0101)
plt.plot(rsq,Pref010101)
plt.plot(rsq,Pref01010101)
plt.legend([r'$\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow$',
r'$\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow$',
r'$\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow$'],fontsize=8)
#plt.title('01010101',fontsize=14)
plt.xlabel("$|r|^2$",fontsize=12)
plt.ylabel("$P_r$",fontsize=12)
plt.grid(b=None, which='major', axis='both')
plt.xlim((0,1))
plt.ylim((0,1.1))
#plt.tight_layout()

plt.subplot(2,2,3)
plt.plot(rsq,Pref1010)
plt.plot(rsq,Pref101010)
plt.plot(rsq,Pref10101010)
plt.legend([r'$\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow$',
r'$\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow$',
r'$\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\uparrow\!\!\!\!\downarrow$'],fontsize=8)
#plt.title('10101010',fontsize=14)
plt.xlabel("$|r|^2$",fontsize=12)
plt.ylabel("$P_r$",fontsize=12)
plt.grid(b=None, which='major', axis='both')
plt.xlim((0,1))
plt.ylim((0,1.1))
#plt.tight_layout()

plt.subplot(2,2,4)
plt.plot(rsq,Pref1100)
plt.plot(rsq,Pref111000)
plt.plot(rsq,Pref11110000)
plt.legend([r'$\uparrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\downarrow$',
r'$\uparrow\!\!\!\!\uparrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\downarrow\!\!\!\!\downarrow$',
r'$\uparrow\!\!\!\!\uparrow\!\!\!\!\uparrow\!\!\!\!\uparrow\!\!\!\!\downarrow\!\!\!\!\downarrow\!\!\!\!\downarrow\!\!\!\!\downarrow$'],fontsize=8)
#plt.title('11110000',fontsize=14)
plt.xlabel("$|r|^2$",fontsize=12)
plt.ylabel("$P_r$",fontsize=12)
plt.grid(b=None, which='major', axis='both')
plt.xlim((0,1))
plt.ylim((0,1.1))
#plt.tight_layout()
ax.set_aspect(1.)
#plt.show()
plt.savefig('PrefPlotSym.jpg',dpi=300)
