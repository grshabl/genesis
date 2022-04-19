# Summary

Our program was written with the main goal as a not intersect gens variants classification between tissues.
The first part of our program describes how to simulate data for testing. 
Second part is about solving our main task.
Third part helps to illustrate our results.

## Main class in our program is Gens.py

p1,p2 - probability of finding sample from the first tissue in the first and in the second moment of corresponding.

| p1           | p2  |
|:-------------| -----:|
| 0.5           | 0.3  |

## Real frequencies
```
Theta [[0.33 0.46 0.37 0.32 0.24 0.4  0.25 0.65 0.77 0.22 0.0.0.0. 0.0.0.0.0.0.]
      [0.0.0.0.0.0.0.0.0.0. 0.53 0.31 0.34 0.7 0.04 0.05 0.02 0.58 0.52 0.62]]
```
## Predicted frequencies 
```
Theta [[0.26240905 0.49145384 0.48545982 0.32293321 0.22608073 0.46093272 0.35151834 0.75688546 0.83183305 0.22423493 0.0.0.0.0. 0.03751189 0.0.0.0.]
      [0.0.0.0.0.0. 0.0.0.0. 0.55203559 0.30672493 0.38596224 0.63180009 0.01666199 0. 0.03342804 0.57037003 0.48352041 0.40781309]]
```

## Example of output:

![plot](https://github.com/grshabl/genesis/blob/main/img/figure1.png)
