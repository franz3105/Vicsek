import matplotlib.pyplot as plt
import numpy as np

rho=[0.8]
for i in range(len(rho)):
  data=np.genfromtxt("output_msd_blocks2.dat", skip_footer=1)

  t=data[:,3]
  MSD=data[:,4]

  plt.plot(t,MSD,label=str(rho[i]))


plt.legend() 
plt.xscale('log')
plt.yscale('log')
plt.xlabel('t')
plt.ylabel('MSD')
plt.title('MSD')
plt.show()