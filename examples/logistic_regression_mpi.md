logistic regression using MPI
=====

This example shows how to use `learn_CRF_mpi` to learn a logistic regression model using MPI.  This uses identical files to the regular [logistic regression](logistic_regression.md) example.

The call is essentially the same as for `learn_CRF`, except prefixed by `mpirun`

```
mpirun ./learn_CRF_mpi -m examples/logistic_regression/model*.txt -d examples/logistic_regression/data*.txt -i 1
```

This eventually results in the output of

```
weights after optimization:
 W[0]:
-0.696163   1.97284
  1.97284 -0.696163
 -1.27668  -1.27668
```