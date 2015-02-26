co_mat
======
covariance matrix calculater.

## example
```ruby
require 'co_mat'

# initialize
# make covariant matrix with 3 samples and 2 dimension
a = CoMat.new([[1,2],[3,4],[5,6]])
a.dim #=> 2
a.sample_count #=> 3

# add sample
a.add([7,8])
a.sample_count #=> 4
a == CoMat.new([[1,2],[3,4],[5,6],[7,8]]) #=> true

# merge sample
b = CoMat.new([[9,9],[11,11]])
a += b
a.sample_count #=> 6

# get covariant matrix
# if samples have integer elements, co_mat return rational value.
a.matrix #=> Matrix[[(35/3), (31/3)], [(31/3), (83/9)]]

# get mean value
a.mean #=> Matrix[[(6/1)], [(20/3)]]

# get autocorrelation matrix
a.r #=> Matrix[[(143/3), (151/3)], [(151/3), (161/3)]]

# linear discriminant analysis solver
# return eigen values and eigen vector matrix
b = CoMat.new([[1,2],[1,3],[2,1]])
lda_solver([b, a]) #=> [Vector[0.7905882352941178, 0.0], Matrix[[-0.354654234120539, 0.7071067811865476], [0.9349975263177833, -0.7071067811865476]]]
```
