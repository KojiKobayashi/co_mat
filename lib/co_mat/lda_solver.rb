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

class LdaSolver

  def self.run(co_mats)
    sw = get_within_matrix(co_mats)
    sb = get_between_matrix(co_mats)
    g_eigen(sb, sw)
  end

  def self.count_all_samples(co_mats)
    co_mats.inject(0) do |count, mat|
      count += mat.sample_count
    end
  end

  def self.get_within_matrix(co_mats)
    all_sample_count = count_all_samples(co_mats)

    co_mats.inject(Matrix.zero(co_mats[0].dim)) do |wMat, co|
      wMat += co.matrix * co.sample_count / all_sample_count
    end
  end

  def self.get_between_matrix(co_mats)
    all_sample_count = count_all_samples(co_mats)

    r = Matrix.zero(co_mats[0].dim)       #autocovariance
    mean = Matrix.zero(co_mats[0].dim, 1) #average

    co_mats.each do |mat|
      r += mat.mean * mat.mean.t * mat.sample_count
      mean += mat.mean * mat.sample_count
    end

    (r - mean * mean.t / all_sample_count) / all_sample_count
  end

  private

  # return eigenvalue vector and eigenvector matrix.
  def self.g_eigen(a, b)
    l = b.cholesky_decomp
    a_dash = l.t.inv * a * l.inv
    eigen = a_dash.eigen

    eigen_set = eigen.eigenvalues.map.each_with_index do |d, i|
      {:eigen_value => d, :eigen_vector => eigen.v.column_vectors[i]}
    end

    # sort by descent of eigen values
    eigen_set.sort!{|a, b| b[:eigen_value].abs <=> a[:eigen_value].abs}

    eigen_value_vec = Vector.elements(eigen_set.map{|e| e[:eigen_value]})
    eigen_vectors = eigen_set.map{|e| e[:eigen_vector]}
    eigen_vectors = eigen_vectors.map{|v| l.inv * v}

    # size normalize
    # TODO: should this normalize be done ?
    eigen_vectors = eigen_vectors.map{|i| i = i / i.norm}

    return eigen_value_vec, Matrix.columns(eigen_vectors)
  end
end
