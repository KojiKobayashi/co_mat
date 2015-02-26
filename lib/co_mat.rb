require "matrix"
sub_dir =File.dirname(__FILE__)
require sub_dir + "/co_mat/lda_solver"

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
