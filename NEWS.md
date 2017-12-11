# symDMatrix 1.0.0.9000

- symDMatrix is now based on the LinkedMatrix package.
- The `symDMatrix` class inherits from `RowLinkedMatrix`.
- Only storing the upper triangular matrix resulted in inefficient queries as
  requests to the lower triangle needed to be redirected. We now store the
  whole matrix, but use virtual transposes for the lower triangular matrix.
  Virtual transposes are very efficient as the block shares the same memory
  mapping as the block across the diagonal and the indices are rewritten
  locally.
- Matrix-like objects that do not support virtual transposes have been dropped
  (i.e., only the ff package is currently left as far as we know).
- `as.symDMatrix()` has been kept the same, but the S4 constructor has changed.


# symDMatrix 1.0.0

Initial release.
