## symDMatrix

### Memory-mapped distributed symmetric matrix

Contact: Gustavo de los Campos (gdeloscampos@gmail.com), Paulino Perez-Rodriguez (perpdgo@gmail.com  )

Class: ```symMatrix``` 

Slots: 
      - names (character)
      - data (list) each element of the list is an ff object
      - centers (numeric) column-means used in the computation of the matrix
      - scales (numeric) column-standard deviations used to scale the matrix.

[Examples]()

Pending:
     - Replacement: function to update content of the matrix
     - Add chunk:   function to add one chunk
     - load2:       a method for loading the object and oppening connecctions
                    (currently load() works provided that you are working on the
                     folder that contains the data). 
                     A template for load2() method can be obtained from BGData