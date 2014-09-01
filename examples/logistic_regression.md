single node inference
=====

This example shows how to use the `infer_CRF` executable to do inference on pretty much the simplest possible graph: a single-node model.

```
   (y0)
```

In this model, we have that

```
p(y0|x) = exp(f0(y0;x))/Z
```

Where `y0` can have the `3` possible values 0,1,2.  The input features are of length 2.  In Marbl, Conditional Random Fields are linearly parameterized.  Thus, if we consider the values

```
f0(y0=0;x)
f0(y0=1;x)
f0(y0=3;x)
```

as a single vector of length 3, these are given by multiplying the input features `f0` with a 3x2 matrix W.

You can run this example by doing, from the command line

```
./learn_CRF -m examples/logistic_regression/model*.txt -f examples/logistic_regression/data*.txt
```

This will read the graph specified in `model.txt` and the parameters specified in `theta.txt`, do inference, and produce marginals in the file `marginals.txt’.  `-i 1` specifies that a single pass of updates over the variables is to be done.


`model1.txt`, `model2.txt`, and `model3.txt` are identical, and have the following form (very similar to that used in the [single node MRF inference](singlenode_inference.md] example.

```
1 3
1 1
1 1
1 0
```

These lines are explained as follows:

1. `1 3` specifies that there is `1` node, that can have `3` possible values.
2. `1 1` specifies that there is `1` factor, with an entropy counting number of `1`
3. `1 0` specifies that there is `1` factor (must match above line) with a "type" of `0`.
4. `1 0` specifies that the first (and only) factor contains `1` node, namely node `0`.

This example will use the following tiny dataset

```
input  1: 0 1
output 1: 0
input  2: 1 0
output 2: 1
input  3 -1 -1
output 3: 2
```

This is specified for Marbl with `data1.txt` having the following form:

```
1
2 0 1
1 0
```

This line is explained as follows:

1. `1`  Specifies that there is `1` factor
2. `2 0 1` Specifies that the first (and only) factor has a feature vector of length `2`, namely `0 1`
3. `1 0` Specifies that there is `1` output node, with value `0`.

`data2.txt` and `data3.txt` are analogous, with the other data.  These are

```
1
2 1 0
1 1
```

```
1
2 -1 -1
1 2
```

When running the code, the following is observed:

```
initial weights:
 W[0]:
0 0
0 0
0 0
calling lbfgs:
Iteration 1:
  fx = 0.586507, x[0] = 0.000000, x[1] = 0.500000
  xnorm = 1.000000, gnorm = 0.371537, step = 1.500000
Iteration 2:
  fx = 0.255431, x[0] = -0.335482, x[1] = 1.341766
  xnorm = 2.418881, gnorm = 0.106625, step = 1.000000
Iteration 3:
  fx = 0.210631, x[0] = -0.548295, x[1] = 1.714221
  xnorm = 3.032677, gnorm = 0.036899, step = 1.000000
Iteration 4:
  fx = 0.203158, x[0] = -0.661671, x[1] = 1.911546
  xnorm = 3.362738, gnorm = 0.008016, step = 1.000000
Iteration 5:
  fx = 0.202745, x[0] = -0.692718, x[1] = 1.966186
  xnorm = 3.454697, gnorm = 0.000861, step = 1.000000
Iteration 6:
  fx = 0.202740, x[0] = -0.696176, x[1] = 1.972667
  xnorm = 3.465690, gnorm = 0.000038, step = 1.000000
Iteration 7:
  fx = 0.202740, x[0] = -0.696163, x[1] = 1.972841
  xnorm = 3.466021, gnorm = 0.000012, step = 1.000000
L-BFGS optimization terminated with status code = 0
  fx = 0.202740, x[0] = -0.696163, x[1] = 1.972841
weights after optimization:
 W[0]:
-0.696163   1.97284
  1.97284 -0.696163
 -1.27668  -1.27668
```

The output weights are written to a file `W.txt` 

```
1
3 2
-0.696163 
1.97284 
-1.27668 
1.97284 
-0.696163 
-1.27668 
```

This is explained as follows:
1. There is only `1` type of factor
2. The first (and only) factor has a weight matrix of size `3` by `2`
3. The values of the matrix are (in column-major order) `-0.696163`, `1.97284`, etc.

Again, this makes intuitive sense.  The second feature is correlated with label `0`, the second feature with label `1`, and both features are negatively correlated with label `2`.