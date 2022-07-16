import numpy as np
import os
import matplotlib.pyplot

test3 = np.ndarray((4,))
print(test3)


test = {'a': 1, 'b': 2, 'c': 3, 'd': 4}

targ_str = 'E1 =     -0.00039027 is the energy of the lowest SO-state'
targ_str = targ_str.replace('E1 =', '')
targ_val = float(targ_str.split(' is')[0])

print(targ_str)
print(targ_val)

blahd = 'hi'
print('u' + blahd)


