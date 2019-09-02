import numpy as np

beta  = np.linspace(1,2,6)
nrand = np.linspace(1,2,6)

counter = 0
for b in beta:
    for n in nrand:
        print(b, n, counter)
        counter +=1

