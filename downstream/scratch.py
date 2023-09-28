import numpy as np




M_ = np.exp(np.array([[-.005,.005],[.95,-.95]]))
for _ in range(100):
    M_ = np.exp(-M_)

import cassiopeia as cs
