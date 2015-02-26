require "matrix"

class Matrix

  # setter
  def []=(i, j, x)
    @rows[i][j] = x
  end

  # cholesky decomposition.
  # return upper triangle matrix.
  def cholesky_decomp
    Matrix.Raise ErrDimensionMismatch unless square?
    a = clone
    dim = a.row_size
    l = Matrix.I(dim)

    (0...dim).each do |i|
      li = Matrix.I(dim)

        li[i,i] = Math::sqrt(a[i,i])
        ((i + 1)...dim).each do |j|
          li[i,j] = 0.0
          li[j,i] = a[j,i] / li[i,i]
        end

        ((i + 1)...dim).each do |m|
          ((i + 1)...dim).each do |n|
            a[m,n] -= a[i,m] * a[n,i] / a[i,i].to_f
          end
        end
        a[i,i] = 1.0
        ((i + 1)...dim).each do |m|
          a[i,m] = 0.0
          a[m,i] = 0.0
        end

      l *= li
    end
    l.t
  end
end

class CoMat

  def initialize(arrs = nil)
    @sample_count = 0
    make_comatrix(arrs) unless arrs == nil
  end

  attr_accessor :sample_count, :xx, :x, :dim

  # initialize and set features.
  # return set samples count.
  def make_comatrix(arrs)
    init_member(arrs[0].count)
    update_mean_and_r(arrs)
  end

  # add one feature and update.
  # return add samples count.
  def add(arr)
    init_member(arr.count) if @sample_count == 0
    update_mean_and_r([arr])
  end

  def ==(other)
    @x == other.x && @xx == other.xx && @sample_count == other.sample_count && @dim == other.dim
  end

  def +(other)
    raise "dimension error #{@dim} != #{other.dim}" if @dim != other.dim
    ret = self.clone
    ret.x += other.x
    ret.xx += other.xx
    ret.sample_count += other.sample_count
    ret
  end

  def matrix
    (r  - mean * mean.t)
  end

  def mean
    @x / @sample_count.to_r
  end

  def r
    @xx / @sample_count.to_r
  end

  private
  def update_mean_and_r(arrs)
    arrs.each do |arr|
      raise "dimension error" unless arr.count == @dim
      vec = Matrix[arr].t
      @x  += vec
      @xx += vec * vec.t
    end
    @sample_count += arrs.count
    arrs.count
  end

  def init_member(dim)
    @dim = dim
    @sample_count = 0
    @x = Matrix.zero(@dim, 1)
    @xx = Matrix.zero(@dim)
  end
end

def lda_solver(co_mats)

  def count_all_samples(co_mats)
    co_mats.inject(0) do |count, mat|
      count += mat.sample_count
    end
  end

  def get_within_matrix(co_mats)
    all_sample_count = count_all_samples(co_mats)

    co_mats.inject(Matrix.zero(co_mats[0].dim)) do |wMat, co|
      wMat += co.matrix * co.sample_count / all_sample_count
    end
  end

  def get_between_matrix(co_mats)
    all_sample_count = count_all_samples(co_mats)

    r = Matrix.zero(co_mats[0].dim)       #autocovariance
    mean = Matrix.zero(co_mats[0].dim, 1) #average

    co_mats.each do |mat|
      r += mat.mean * mat.mean.t * mat.sample_count
      mean += mat.mean * mat.sample_count
    end

    (r - mean * mean.t / all_sample_count) / all_sample_count
  end

  # return eigenvalue vector and eigenvector matrix.
  def g_eigen(a, b)
    l = b.cholesky_decomp
    a_dash = l.t.inv * a * l.inv
    eigen = a_dash.eigen

    eigen_set = eigen.eigenvalues.map.each_with_index do |d, i|
      {:eigen_value => d, :eigen_vector => eigen.v.column_vectors[i]}
    end

    # sort by descent of eigen values
    eigen_set.sort!{|a, b| b[:eigen_value] <=> a[:eigen_value]}

    eigen_value_vec = Vector.elements(eigen_set.map{|e| e[:eigen_value]})
    eigen_vectors = eigen_set.map{|e| e[:eigen_vector]}
    eigen_vectors = eigen_vectors.map{|v| l.inv * v}

    # size normalize
    # TODO: should this normalize be done ?
    eigen_vectors = eigen_vectors.map{|i| i = i / i.norm}

    return eigen_value_vec, Matrix.columns(eigen_vectors)
  end

  sw = get_within_matrix(co_mats)
  sb = get_between_matrix(co_mats)
  g_eigen(sb, sw)
end

