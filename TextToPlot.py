import numpy as np
from matplotlib import pyplot as plt

Pref4 = np.loadtxt("RTnucspin4.txt")
Pref6 = np.loadtxt("RTnucspin6.txt")
Pref8 = np.loadtxt("RTnucspin8.txt")
Pref10 = np.loadtxt("RTnucspin10.txt")

r = np.linspace(0,1,101)
rsq = r**2

plt.figure(1)
plt.subplots_adjust(left=.1, bottom=.12, right=.97, top=.94)

plt.plot(rsq,Pref4,label='4 Nuclei')
plt.plot(rsq,Pref6,label='6 Nuclei')
plt.plot(rsq,Pref8,label='8 Nuclei')
plt.plot(rsq,Pref10,label='10 Nuclei')

plt.legend()
ax = plt.gca()
plt.legend(['4 Nuclei','6 Nuclei','8 Nuclei','10 Nuclei'], bbox_to_anchor=(1,.325), bbox_transform=ax.transAxes, fontsize=12)
plt.title('Probability of Total Reflection')
plt.xlabel("$|r|^2$",fontsize=14)
plt.ylabel("$P_r$",fontsize=14)
plt.grid(b=None, which='major', axis='both')
plt.xlim((0,1))
plt.ylim((0,.9))
#plt.show()
plt.savefig('PrefPlots.jpg',dpi=300)