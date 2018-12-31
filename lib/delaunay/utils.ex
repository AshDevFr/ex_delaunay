defmodule Delaunay.Utils do
  use Bitwise

  @moduledoc """
  Documentation for Delaunay.Utils
  """

  @doc """
  pseudoAngle: Monotonically increases with real angle, but doesn't need expensive trigonometry
  """
  def pseudoAngle(dx, dy) do
    p = dx / (abs(dx) + abs(dy))
    (if (dy > 0), do: 3 - p, else: 1 + p) / 4 # [0..1]
  end

  @doc """

  """
  def dist(ax, ay, bx, by) do
    dx = ax - bx
    dy = ay - by
    dx * dx + dy * dy
  end

  @doc """

  """
  def orient(px, py, qx, qy, rx, ry) do
    (qy - py) * (rx - qx) - (qx - px) * (ry - qy) < 0
  end

  @doc """

  """
  def inCircle(ax, ay, bx, by, cx, cy, px, py) do
    dx = ax - px
    dy = ay - py
    ex = bx - px
    ey = by - py
    fx = cx - px
    fy = cy - py

    ap = dx * dx + dy * dy
    bp = ex * ex + ey * ey
    cp = fx * fx + fy * fy

    dx * (ey * cp - bp * fy) -
      dy * (ex * cp - bp * fx) +
      ap * (ex * fy - ey * fx) < 0
  end

  @doc """

  """
  def circumradius(ax, ay, bx, by, cx, cy) do
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    d = if ((dx * ey - dy * ex) == 0), do: :math.pow(2, 50), else: 0.5 / (dx * ey - dy * ex)

    x = (ey * bl - dy * cl) * d
    y = (dx * cl - ex * bl) * d

    x * x + y * y
  end

  @doc """

  """
  def circumcenter(ax, ay, bx, by, cx, cy) do
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    d = if ((dx * ey - dy * ex) == 0), do: :math.pow(2, 50), else: 0.5 / (dx * ey - dy * ex)

    x = ax + (ey * bl - dy * cl) * d
    y = ay + (dx * cl - ex * bl) * d

    {x, y}
  end

  @doc """

  """
  def quicksort(ids, dists, left, right) do
    if (right - left <= 20) do
      List.foldl(
        ((left + 1)..right)
        |> Enum.to_list,
        ids,
        fn i, ids ->
          tmp = ids
                |> Enum.at(i)
          tmp_dist = dists
                     |> Enum.at(tmp)

          swipeRight(ids, dists, tmp_dist, left, i - 1)
          |> (
               fn {j, l} ->
                 l
                 |> List.replace_at(j + 1, tmp)
               end).()
        end
      )
    else
      median = (left + right)
               |> bsr(1)

      i = left + 1

      new_ids = ids
                |> swap(median, i)
                |> (fn l -> if (dist_x(l, dists, left) > dist_x(l, dists, right)) do
                              l
                              |> swap(left, right)
                            else
                              l
                            end
                    end).()
                |> (fn l -> if (dist_x(l, dists, i) > dist_x(l, dists, right)) do
                              l
                              |> swap(i, right)
                            else
                              l
                            end
                    end).()
                |> (fn l -> if (dist_x(l, dists, left) > dist_x(l, dists, i)) do
                              l
                              |> swap(left, i)
                            else
                              l
                            end
                    end).()

      tmp = new_ids
            |> Enum.at(i)
      tmp_dist = dists
                 |> Enum.at(tmp)

      new_ids
      |> swapIds(dists, tmp_dist, i + 1, right - 1)
      |> (fn {i, j, l} ->
        {
          i,
          j,
          l
          |> List.replace_at(
               left + 1,
               l
               |> Enum.at(j)
             )
          |> List.replace_at(j, tmp)
        }
          end).()
      |> (fn {i, j, l} ->
        if (right - i + 1 >= j - left) do
          l
          |> quicksort(dists, i, right)
          |> quicksort(dists, left, j - 1)
        else
          l
          |> quicksort(dists, left, j - 1)
          |> quicksort(dists, i, right)
        end
          end).()
    end
  end

  @doc """

  """
  def dist_x(ids, dists, x) do
    dists
    |> Enum.at(
         ids
         |> Enum.at(x)
       )
  end

  @doc """

  """
  def swapIds(ids, dists, tmp_dist, i, j) do
    cond do
      dist_x(ids, dists, i) < tmp_dist ->
        ids
        |> swapIds(dists, tmp_dist, i + 1, j)
      dist_x(ids, dists, j) > tmp_dist ->
        ids
        |> swapIds(dists, tmp_dist, i, j - 1)
      j < i ->
        {i, j, ids}
      true ->
        ids
        |> swap(i, j)
        |> swapIds(dists, tmp_dist, i, j)

    end
  end

  @doc """

  """
  def swipeRight(ids, dists, tmp_dist, left, j) do
    if (j >= left && Enum.at(dists, Enum.at(ids, j)) > tmp_dist) do
      swipeRight(
        ids
        |> List.replace_at(
             j + 1,
             ids
             |> Enum.at(j)
           ),
        dists,
        tmp_dist,
        left,
        j - 1
      )
    else
      {j, ids}
    end
  end

  @doc """

  """
  def swap(list, i, j) do
    tmp = list
          |> Enum.at(i)
    list
    |> List.replace_at(
         i,
         (
           list
           |> Enum.at(j))
       )
    |> List.replace_at(j, tmp)
  end

  @doc """

  """
  def defaultGetX({x, _y}) do
    x
  end

  @doc """

  """
  def defaultGetY({_x, y}) do
    y
  end
end
