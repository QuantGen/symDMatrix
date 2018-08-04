## Test environments

* Local Arch Linux install: R 3.5.1
* Travis CI: R 3.5.0
* win-builder: R Under development (unstable)


## R CMD check results

There were 0 errors, 1 warning, and 0 notes.

### WARNING: '::' or ':::' import not declared from: 'BEDMatrix'

This package depends on the BGData package which depends on the BEDMatrix
package. Therefore, the import should not need to be declared.


## revdep_check results

There were 0 packages with problems.
