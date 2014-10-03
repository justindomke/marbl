marbl Frequently Asked Question(s)
=====

# Hidden Variables

  * Q) What label should I use for a node whose value is hidden?
  * A) Valid labels start at 0.  Thus, hidden labels are declared simply by using a value of -1 for any components that are not observed.  When computing the marginal-based loss, these nodes will not be used.  (The input features are still used for inference, of course.)

# Contact

  * Please email the author using the contact info [here](http://users.cecs.anu.edu.au/~jdomke/) with any further questions.