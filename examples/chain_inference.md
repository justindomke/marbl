marbl
=====

This example shows how to use the `infer_MRF` executable to do inference on a simple chain model.

```
(x₀)—(x₁)—(x₂)—(x₃)
```

In this problem, we will use the following distribution


```
p(x)=exp(f₀(x₀)+f₁(x₁)+f₂(x₂)+f₃(x₃)+f₀₁(x₀,x₁)+f₁₂(x₁,x₂)+f₁₂(x₂,x₃))/Z
f₀(x₀=0)=-1
f₀(x₀=1)=+1
f₁(x₁=0)=-1
f₁(x₁=1)=+1
f₂(x₂=0)=-1
f₂(x₂=1)=+1
f₃(x₃=0)=-1
f₃(x₃=1)=+1
f₀₁(x₀≠x₁)=-1
f₀₁(x₀=x₁)=+1
f₁₂(x₁≠x₂)=-1
f₁₂(x₁=x₂)=+1
f₂₃(x₂≠x₃)=-1
f₂₃(x₂=x₃)=+1
```

for some normalization factor Z.  You can run this example by doing, from the command line

```
./infer_MRF -m examples/chain_inference/model.txt -f examples/chain_inference/theta.txt -mu marginals.txt -i 1
```

This will read the graph specified in `model.txt` and the parameters specified in `theta.txt`, do inference, and produce marginals in the file `marginals.txt`.  `-i 1` specifies that a single pass of updates over the variables is to be done.  (A single pass updates all cliques from start to end and then in reverse.)

`model.txt` has the following form:

```
4 2 2 2 2
7 0 1 -1 1 -1 1 0
7 0 0 0 0 0 0 0
1 0
2 0 1
1 1
2 1 2
1 2
2 2 3
1 3
```

These lines are explained as follows:

1. `4 2 2 2 2` specifies that there are `4` nodes, each of which can have `2` possible values.
2. `7 0 1 -1 1 -1 1 0 ` specifies that there are `7` factors, with an entropy counting numbers of `0`, `1`, etc.  (More on this below)
3. `7 0 0 0 0 0 0 0`  specifies that there are `7` factors (must match above line) with a "type" of `0`.  (This is not used with MRF inference, but must be specified to be consistent with the model format used with CRFs)
4. `1 0` specifies that the first factor contains `1` node, namely node `0`.
5. `2 0 1` specifies that the second factor contains `2` nodes, namely `0` and `1`.
6. `1 1` specifies that the third factor contains `1` node, namely node `1`.
7. `2 1 2` specifies that the fourth factor contains `2` nodes, namely `1` and `2`.
8. `1 2` specifies that the fifth factor contains `1` node, namely node `2`.
9. `2 2 3` specifies that the sixth factor contains `2` nodes, namely `2` and `3`.
10. `1 3` specifies that the seventh factor contains `1` node, namely node `3`.

Why put the factors in this order?  This is because inference in this toolbox proceeds in the order specified here.  Namely, one iterations corresponds to, first, updating all factors in the order given, and then updating all factors in the reverse of the order given.  For a chain, this means that all messages will have converged in a single iteration

Where do the entropy factors come from?  These correspond to the [Bethe entropy](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.21.2420).  For pairwise graphs, the Bethe entropy gives a count of 1 to all pairs, and a count of 1-(#pairs-touching) for all nodes.  This is exact in treelike graphs, such as above.  However, with high-treewidth graphs, one will only produce approximate marginals, and the selection of entropy constants is still a subject of active research.  Thus, the toolbox gives you the power to select these yourself.

`theta.txt` has the following form:

```
2 -1 1
4 1 -1 -1 1
2 -1 1
4 1 -1 -1 1
2 -1 1
4 1 -1 -1 1
2 -1 1
```

This line is explained as follows:

1. `2 -1 1`  Specifies that the first factor (that for x₀ alone) has `2` possible configurations, and that these have parameter values `-1` and `+1`.
2. `4 1 -1 -1 1` Specifies that the second factor (that for x₀ and x₁) has `4` possible configurations, and that these have parameter values of `1` when x₀=x₁ and `-1` when x0≠x1.
3. The other lines are similar.

Note that the factors for pairs are given in column-major order.  Similarly, if one has a factor with N variables, this is reduced to a single line by going through the first dimension most quickly, and the last dimension most slowly.

After running the code, a file `marginals.txt` will be produced in the main directory:

```
2 0.0209241 0.979076 
4 0.00304756 0.00304756 0.0178766 0.976028 
2 0.00609512 0.993905 
4 0.0012799 0.00481522 0.00481522 0.98909 
2 0.00609512 0.993905 
4 0.00304756 0.0178766 0.00304756 0.976028 
2 0.0209241 0.979076
```

This is explained as follows:
1. The first factor (that for x₀ alone) has `2` possible configurations, and these have a marginal probability of `0.0209241` and `0.979076`
2. The second factor (that for x₀ and x₁) has `4` possible configurations, and these have a marginal probability of `0.00304756 0.00304756 0.0178766 0.976028`.
3. The other lines are similar.

As above, this is in column-major format.  So, for example

```
p(x₀=0,x₁=0)=0.00304756
p(x₀=1,x₁=0)=0.00304756
p(x₀=0,x₁=1)=0.0178766
p(x₀=1,x₁=1)=0.976028 
```

Intuitively, the marginal probabilities are mostly concentrated around an all-one configuration.  This makes sense, since all the unary factors point in that direction, and there is a strong smoothness penalty in the pairwise terms.

You can verify, by running the same code with, e.g. `-i 2` that the marginals do not change, showing that inference has completed after a single pass.