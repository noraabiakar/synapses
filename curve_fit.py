import ampa_cell as ac

from numpy import array, exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# observations
gt = ac.run_ground_truth()
y = gt[:,1]
x = array(range(len(y)))

def pred(x, a, b, c):
    approx = ac.run_approximation(a, b, c)
    return approx[:,1]

params, _ = curve_fit(pred, x, y, bounds=[[-10, -10, -10], [10, 10, 10]])
print(params)
