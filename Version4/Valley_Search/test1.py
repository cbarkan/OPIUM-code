import numpy as np
a = np.array([1,2,3])
f = open('test_text.txt','w')
f.write(np.array2string(a))
f.close()