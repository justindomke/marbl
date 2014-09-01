marbl
=====

This example shows how to use the `infer_MRF` executable to do inference on pretty much the simplest possible graph: a single-node model.  You can run this example by doing, from the command line

```
./infer_MRF -m examples/singlenode_inference/model.txt -f examples/singlenode_inference/theta.txt -mu marginals.txt
```

This will read the graph specified in `model.txt` and the parameters specified in `theta.txt`, do inference, and produce marginals in the file `marginals.txt'

model.txt has the following form:

```
1 2
1 1
1 0
1 0
```

These lines are explained as follows:

First line: `1 2` specifies that there is `1` node, that can have `2` possible values.
Second line: `1 1` specifies that there is `1` factor, with an entropy countring number of `1`
Third line: `1 0`  specifies that there is `1` factor (must match above line) with a "type" of `0`.  (This is not used with MRF inference, but must be specified to be consistent with the model format used with CRFs)
Fourth line: `1 0` specifies that the first (and only) factor contains `1` node, namely node `0`.