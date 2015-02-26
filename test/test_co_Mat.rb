require "test/unit"
require "../lib/co_mat.rb"

class TC_COMAT < Test::Unit::TestCase
  def setup
    @a = [[1,0],[3,2],[2,-1],[-1,-2],[0,1],[9,3]]
    @b = [[-1,0],[1,2],[0,-1],[-3,-2],[-2,1]]
    @ex = [[-1,1,2],[-2,3,1],[4,0,3]]
    @ex_co_matrix = CoMat.new(@ex)
  end

  def test_matrix

    # factor for define of div by (sample count) or (sample count - 1)
    factor = (@ex.length.to_f - 1) / @ex.length
    assert_in_delta(10.3333 * factor, @ex_co_matrix.matrix[0,0], 0.001)
    assert_in_delta(-4.1667 * factor, @ex_co_matrix.matrix[0,1], 0.001)
    assert_in_delta(3.00000 * factor, @ex_co_matrix.matrix[0,2], 0.001)
    assert_in_delta(-4.1667 * factor, @ex_co_matrix.matrix[1,0], 0.001)
    assert_in_delta(2.33333 * factor, @ex_co_matrix.matrix[1,1], 0.001)
    assert_in_delta(-1.5000 * factor, @ex_co_matrix.matrix[1,2], 0.001)
    assert_in_delta(3.00000 * factor, @ex_co_matrix.matrix[2,0], 0.001)
    assert_in_delta(-1.5000 * factor, @ex_co_matrix.matrix[2,1], 0.001)
    assert_in_delta(1.00000 * factor, @ex_co_matrix.matrix[2,2], 0.001)
  end

  def test_equal
    set1 = CoMat.new(@a)
    set2 = CoMat.new(@a)
    assert_equal(set1, set2)

    set2.add([1,2])
    assert_not_equal(set1, set2)

    set1.add([1,2])
    assert_equal(set1, set2)
  end

  def test_merge
    set1 = CoMat.new(@a + @b)
    set2 = CoMat.new(@a)
    set3 = CoMat.new(@b)

    set2 += set3
    assert_equal(set1, set2)
  end

  def test_LDA
    set1 = CoMat.new([[1,0],[3,2],[2,-1],[-1,-2],[0,1]])
    set2 = CoMat.new([[-1,0],[1,2],[0,-1],[-3,-2],[-2,1]])

    d, v = LdaSolver.run([set1, set2])

    gt_v1 = Vector[5.0/Math.sqrt(34), -3.0/Math.sqrt(34)]
    gt_v2 = Vector[0.0, 1.0]

    assert_in_delta(25.0/32, d[0], 0.0001)
    assert_in_delta(0,       d[1], 0.0001)

    # (v = gt) || (v = -gt)
    # => (v // gt) && (v.abs == gt.abs == 1)
    assert_in_delta(1.0, gt_v1.inner_product(v.column_vectors[0]).abs, 0.0001)
    assert_in_delta(1.0, gt_v2.inner_product(v.column_vectors[1]).abs, 0.0001)
    assert_in_delta(1.0, v.column_vectors[0].inner_product(v.column_vectors[0]).abs, 0.0001)
    assert_in_delta(1.0, v.column_vectors[1].inner_product(v.column_vectors[1]).abs, 0.0001)

    # when 2 class,
    # 1st eigen value > 0
    # 2nd eigen value == 0
    set1 = CoMat.new([[1,2],[3,4],[5,6],[7,8]])
    set2 = CoMat.new([[1,2],[2,1],[1,2]])
    d,v = LdaSolver.run([set1, set2])
    assert_equal(true, d[0] > 0.0)
    assert_in_delta(0, d[1], 0.0001)
  end

end
