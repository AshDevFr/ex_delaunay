defmodule Delaunay.Internals do
  alias Delaunay
  alias Delaunay.Utils

  @moduledoc """
  Documentation for Delaunay.Internals
  """

  @doc """

  """
  def hashKey(%Delaunay{hash_size: hash_size, cx: cx, cy: cy}, x, y) do
    Utils.pseudoAngle(x - cx, y - cy) * hash_size
    |> trunc
    |> rem(hash_size)
  end

  @doc """

  """
  def legalize(
        %Delaunay{
          triangles: triangles,
          coords: coords,
          halfedges: halfedges,
          edge_stack: edge_stack,
          hull_start: hull_start
        } = delaunay,
        a,
        i \\ 0
      ) do
    # if the pair of triangles doesn't satisfy the Delaunay condition
    # (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
    # then do the same check/flip recursively for the new pair of triangles
    #
    #           pl                    pl
    #          /||\                  /  \
    #       al/ || \bl            al/    \a
    #        /  ||  \              /      \
    #       /  a||b  \    flip    /___ar___\
    #     p0\   ||   /p1   =>   p0\---bl---/p1
    #        \  ||  /              \      /
    #       ar\ || /br             b\    /br
    #          \||/                  \  /
    #           pr                    pr

    b = halfedges
        |> Enum.at(a)
    a0 = a - rem(a, 3)
    ar = a0 + rem(a + 2, 3)

    if (b == -1) do
      if (i == 0) do
        {delaunay, ar}
      else
        legalize(
          delaunay,
          edge_stack
          |> Enum.at(i - 1),
          i - 1
        )
      end
    else
      b0 = b - rem(b, 3)
      al = a0 + rem(a + 1, 3)
      bl = b0 + rem(b + 2, 3)

      p0 = triangles
           |> Enum.at(ar)
      pr = triangles
           |> Enum.at(a)
      pl = triangles
           |> Enum.at(al)
      p1 = triangles
           |> Enum.at(bl)

      illegal = Utils.inCircle(
        coords
        |> Enum.at(2 * p0),
        coords
        |> Enum.at(2 * p0 + 1),
        coords
        |> Enum.at(2 * pr),
        coords
        |> Enum.at(2 * pr + 1),
        coords
        |> Enum.at(2 * pl),
        coords
        |> Enum.at(2 * pl + 1),
        coords
        |> Enum.at(2 * p1),
        coords
        |> Enum.at(2 * p1 + 1)
      )

      if (illegal) do
        new_triangles =
          triangles
          |> List.replace_at(a, p1)
          |> List.replace_at(b, p0)

        hbl = halfedges
              |> Enum.at(bl)

        # edge swapped on the other side of the hull (rare); fix the halfedge reference
        new_delaunay =
          %{delaunay | triangles: new_triangles}
          |> (fn d ->
            if (hbl == -1) do
              hull_tri(delaunay, hull_start, a, bl)
            else
              d
            end
              end).()
          |> link(a, hbl)
          |> link(
               b,
               halfedges
               |> Enum.at(ar)
             )
          |> link(ar, bl)

        br = b0 + rem(b + 1, 3)

        # don't worry about hitting the cap: it can only happen on extremely degenerate input
        new_delaunay
        |> (fn d ->
          if (
               i < (
                 edge_stack
                 |> length)) do
            {
              i + 1,
              %{
                d |
                edge_stack:
                  edge_stack
                  |> List.replace_at(i, br)
              }
            }
          else
            {i, d}
          end
            end).()
        |> (fn {i, d} ->
          legalize(
            d,
            a,
            i
          )
            end).()
      else
        if (i == 0) do
          {delaunay, ar}
        else
          legalize(
            delaunay,
            edge_stack
            |> Enum.at(i - 1),
            i - 1
          )
        end
      end
    end
  end

  @doc """

  """
  def hull_tri(%Delaunay{hull_start: hull_start, hull_tri: hull_tri, hull_next: hull_next} = delaunay, e, a, bl) do
    if (
         hull_tri
         |> Enum.at(e) == bl) do
      %{
        delaunay |
        hull_tri:
          hull_tri
          |> List.replace_at(e, a)
      }
    else
      new_e =
        hull_next
        |> Enum.at(e)
      if (new_e != hull_start) do
        hull_tri(delaunay, new_e, a, bl)
      else
        delaunay
      end
    end
  end

  @doc """

  """
  def link(%Delaunay{halfedges: halfedges} = delaunay, a, b) do
    new_halfedges = halfedges
                    |> List.replace_at(a, b)
                    |> (fn l -> if (b != -1) do
                                  l
                                  |> List.replace_at(b, a)
                                else
                                  l
                                end
                        end).()
    %{delaunay | halfedges: new_halfedges}
  end

  @doc """
  addTriangle: Add a new triangle given vertex indices and adjacent half-edge ids
  """
  def addTriangle(%Delaunay{triangles_len: triangles_len, triangles: triangles} = delaunay, i0, i1, i2, a, b, c) do
    t = triangles_len
    new_triangles = triangles
                    |> List.replace_at(triangles_len, i0)
                    |> List.replace_at(triangles_len + 1, i1)
                    |> List.replace_at(triangles_len + 2, i2)

    new_triangles_len = triangles_len + 3


    new_delaunay =
      %{delaunay | triangles: new_triangles, triangles_len: new_triangles_len}
      |> link(triangles_len, a)
      |> link(triangles_len + 1, b)
      |> link(triangles_len + 2, c)
    {new_delaunay, t}
  end
end
