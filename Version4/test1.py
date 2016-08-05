import numpy as np
import scipy as sp
test_A = np.array([[ 3.53738667, -1.83163503, -1.04345406],[-1.83163503 , 7.11373881 , 0.1887602 ],[-1.04345406,  0.1887602 ,  1.79958654]])
test_B = np.array([4,-1,3.1])
x = sp.sparse.linalg.cg(test_A,test_B)
print x