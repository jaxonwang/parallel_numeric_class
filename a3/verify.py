import numpy as np

raw_a = [
    [9.68346, 0],
    [0.840188, 10.1806],
]
in_str = \
"""
5.00471 0.840188
0.840188 5.70498
"""

in_str1 = \
"""
0.394383 0.783099
0.79844 0.911647
"""

def str_to_array(in_str):
    raw_b = [[float(j) for j in i.split(" ")] for i in in_str.split('\n')[1:-1]]
    return raw_b;

A = np.array(str_to_array(in_str))
B = np.array(str_to_array(in_str1))
L = np.linalg.cholesky(A)
print(np.linalg.inv(L.T))
print(np.matmul(B, np.linalg.inv(L.T)))

in_str =\
"""
5.00471 0 0 0 0 0 0 0
0.840188 5.70498 0 0 0 0 0 0
0.394383 0.783099 4.53907 0 0 0 0 0
0.79844 0.911647 0.197551 4.27329 0 0 0 0
0.335223 0.76823 0.277775 0.55397 4.89534 0 0 0
0.477397 0.628871 0.364784 0.513401 0.95223 4.08277 0 0
0.916195 0.635712 0.717297 0.141603 0.606969 0.0163006 4.14288 0
0.242887 0.137232 0.804177 0.156679 0.400944 0.12979 0.108809 2.98052
"""
Lall = np.linalg.cholesky(np.array(str_to_array(in_str)))

print(Lall)
