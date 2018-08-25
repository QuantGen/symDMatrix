## Test environments

* Local Arch Linux install: R 3.5.1
* Travis CI: R 3.5.0
* win-builder: R Under development (unstable)


## R CMD check results

There were 0 errors, 0 warning, and 0 notes.


## revdep_check results

There were 0 packages with problems.


## Comments

Sorry for the frequent updates. I'm trying to fix the CRAN check issue by
introducing a read-only mode to the package and figuring out the places where I
need to enable it to pass the check. Note that I am using `ff`, which I don't
think actually writes to the user library, but probably opens the example data
files for writing.
