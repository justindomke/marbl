marbl
=====

This example shows how to use the `infer_MRF` executable to do inference on a simple chain model.  You can run this example by doing, from the command line

  ./infer_MRF -m examples/chain_inference/model.txt -f examples_chain_inference/theta.txt -mu marginals.txt

This will read the graph specified in `model.txt` and the parameters specified in `theta.txt`, do inference, and produce marginals in the file `marginals.txt'