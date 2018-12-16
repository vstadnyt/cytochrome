60-Hz-noise
=====

Summary:
-----
A small tool for calculating 60-Hz in time-resolved data

Authors:
-----
Valentyn Stadnytskyi

Available background functions:
-----
* fixed freqency noise subtraction

Example Usage:
-----
```python
import numpy as np
from relax import relaxation_fit, single_step_relaxation

x = [1,2,5,10,15,25,35,60,90, 200, 500, 1000, 10000, 1000000, 10000000000]
y = [25*(1-np.exp(-5*i))+2 for i in x]

parameters, covariances, y_calc = relaxation_fit(x, y, relaxation_function = single_step_relaxation, initial_guess=[18, 11, 10])

print(parameters) 
# Ideally this should converge to 25., 5., 2. for this example - more data points will improve convergence.

```

Requirements:
-----
* Python >= 3.6
* Numpy
* Scipy
