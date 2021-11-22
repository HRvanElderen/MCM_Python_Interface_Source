import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MSC

f = open("MinCompSpin-main/INPUT/Dataset_Shapes_n9_N1e5.dat", 'r')

data = [[],[],[],[],[],[],[],[],[]]
for line in f:
    for i, v in enumerate(list(line)):
        if v == '\n':
            continue
        data[i].append(int(v))

np_data = np.array(data)
corr = np.corrcoef(np_data)
#data = pd.DataFrame(data) 
#corr = data.corr()
print(corr)
fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(corr,cmap='coolwarm', vmin=-1, vmax=1)
fig.colorbar(cax)
ticks = np.arange(0,9,1)
ax.set_xticks(ticks)
plt.xticks(rotation=90)
ax.set_yticks(ticks)
plt.show()