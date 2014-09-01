single node inference
=====

This example shows how to use the `infer_MRF` executable to do inference on pretty much the simplest possible graph: a single-node model.

```
   (x0)
```

In this model, we have that

```
p(x0=0) = exp(-1)/Z
p(x0=1) = exp(+1)/Z
```

for some normalization factor Z.  You can run this example by doing, from the command line

```
./infer_MRF -m examples/singlenode_inference/model.txt -f examples/singlenode_inference/theta.txt -mu marginals.txt -i 1 
```

This will read the graph specified in `model.txt` and the parameters specified in `theta.txt`, do inference, and produce marginals in the file `marginals.txt’.  `-i 1` specifies that a single pass of updates over the variables is to be done.


`model.txt` has the following form:

```
1 2
1 1
1 0
1 0
```

These lines are explained as follows:

1. `1 2` specifies that there is `1` node, that can have `2` possible values.
2. `1 1` specifies that there is `1` factor, with an entropy countring number of `1`
3. `1 0`  specifies that there is `1` factor (must match above line) with a "type" of `0`.  (This is not used with MRF inference, but must be specified to be consistent with the model format used with CRFs)
4. `1 0` specifies that the first (and only) factor contains `1` node, namely node `0`.

`theta.txt` has the following form:

```
2 -1 1
```

This line is explained as follows:

1. `2 -1 1`  Specifies that the first (and only) factor has `2` possible configurations, and that these have parameter values `-1` and `+1`.

After running the code, a file `marginals.txt` will be produced in the main directory:

```
2 0.119203 0.880797 
```

This is explained as follows:
1. The first (and only) factor has `2` possible configurations, and these have a marginal probability of `0.119203` and `0.880797`

Note that this makes sense.  If we look at the original specification of the model, we could have calculated the probabilities by hand:

```
p(x0=0)=exp(-1)/(exp(-1)+exp(+1))=0.119203
p(x1=1)=exp(+1)/(exp(-1)+exp(+1))=0.880797
```