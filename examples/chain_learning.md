Chain CRF learning
=====

This example shows how to use the `learn_CRF` executable to learn CRF weights for a simple chain model with tied parameters

```
   (y0)—(y1)—(y2)—(y3)
```

This example uses a distribution similar to that in the [chain inference](chain_inference.md) example.  You can run this example by doing, from the command line

```
./learn_CRF -m examples/chain_learning/model*.txt -d examples/chain_learning/data*.txt
```

In this example, there are two types of factors: those for single-nodes, and those for pairs.  This is specified in `model1.txt` (`model2.txt` and `model3.txt` are identical).

```
4 2 2 2 2
7 0 1 -1 1 -1 1 0
7 0 1  0 1  0 1 0
1 0
2 0 1
1 1
2 1 2
1 2
2 2 3
1 3
```

All lines other than the third are exactly the same as in the [chain inference](chain_inference.md) example.  The only difference is the third line

3. `7 0 1  0 1  0 1 0` states that there are 7 factors, and these have types alternating between `0` and `1`.  If you look at the subsequent lines, you can see that the type `0` factors correspond to single nodes, while the type `1` factors correspond to pairs.

The data files specify the input features.  For example, `data1.txt` contains

```
7
2 -.1 .9
1 .9
2 .1 .7
1 .7
2 .1 .6
1 0.1 
2 .1 -.5
4 1 1 1 0
```

The lines are explained by stating that there are `7` factors. The first factor has `2` input features, namely, `-.1` and `.9`.  The second factor has `1` input feature, namely `.9`.  Finally, the last line states that there are `4` output nodes, with labels `1 1 1 0`, respectively.

`data2.txt` and `data3.txt` are similar:

```
7
2 -.1 -.5
1 .05
2 .1 .7
1 .7
2 -.1 .8
1 0.9 
2 .1 -.5
4 0 1 1 1
```

```
7
2 -.1 .9
1 .9
2 .1 .7
1 .1
2 -.1 .01
1 0.9 
2 .1 -.5
4 1 1 0 0
```

If you look at these datasets closely, you can observe the following features:

1. The first input feature for the single-node factors is essentially random.
2. The second input feature is correlated with the corresponding output node being `1`, and negatively correlated with it being `0`.
3. The (single) input feature for pair nodes is correlated with the variables in that node being the same.

If we run the code, learning eventually produces the output of

```
 W[0]:
-0.174765  -1.66048
 0.174765   1.66048
 W[1]:
 0.899209
 -1.51373
-0.312409
 0.926929
```

This makes intuitive sense, given the correlations discussed above:

1. `W[0]` strongly encourages the corresponding output label to be `1` then the second input feature is positive, and largely ignores the first input
2. `W[1]` encourages the labels to be the same when the input feature is large.

This is written to a file `W.txt`.

```
2
2 2
-0.174765 
0.174765 
-1.66048 
1.66048 
4 1
0.899209 
-1.51373 
-0.312409 
0.926929 
```

This is easy to understand.

1. The first line states that there are `2` weight matrices.
2. The second line states that the first matrix is of size `2` by `2`.
3. The third to sixth lines gives the entries of the first matrix
4. The seventh line states that the second matrix is of size `4` by `1`.
8. The tenth to eleventh lines give the entries of the second matrix.

# Inference

Once the model has been learned, one might wish to use it to predict marginals on some test datum.  This can be done by running

```
./infer_CRF -m examples/chain_learning/model1.txt -d examples/chain_learning/testdata.txt -w W.txt -mu marginals.txt -i 1
```

The file `testdata.txt` contains

```
7
2 -.1 .9
1 .1
2 .1 -.1
1 1.3
2 -.1 -.8
1 0.2 
2 .1 1.1
```

Intuitively, we would expect nodes `0` and `3` to be pushed towards output `1`, and nodes `1` and `2` to be pushed towards output `0`.  After running inference, we can examine the file `marginals.txt`:

```
2 0.0602732 0.939727 
4 0.0555329 0.836806 0.00474038 0.102921 
2 0.892339 0.107661 
4 0.876833 0.0354171 0.0155058 0.0722443 
2 0.91225 0.0877501 
4 0.0281865 0.00132725 0.884063 0.0864228 
2 0.0295138 0.970486 
```

This has the factors in the order specified above.  We can see that the marginals probabilities fit exactly what we expect.

