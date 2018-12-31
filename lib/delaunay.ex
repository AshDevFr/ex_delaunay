defmodule Delaunay do
  use Bitwise
  alias __MODULE__
  alias Delaunay.Internals
  alias Delaunay.Utils

  @moduledoc """
  Documentation for Delaunay.
  """

  @epsilon :math.pow(2, -52)

  defstruct [
    :coords,
    :halfedges,
    :hull,
    :hull_next,
    :hull_prev,
    :hull_start,
    :hull_tri,
    :hull_hash,
    :hull_size,
    :triangles,
    :triangles_len,
    :hash_size,
    :cx,
    :cy,
    :edge_stack
  ]

  @doc """

  """
  def from(points) do
    n = length(points)
    {_, coords} =
      points
      |> List.foldl(
           {0, List.duplicate(0, n * 2)},
           fn p, {i, c} ->
             {
               i + 1,
               c
               |> List.replace_at(
                    2 * i,
                    p
                    |> Utils.defaultGetX
                  )
               |> List.replace_at(
                    2 * i + 1,
                    p
                    |> Utils.defaultGetY
                  )
             }
           end
         )

    new(coords)
  end

  @doc """

  """
  def new(coords) do
    n = length(coords)
        |> bsr(1)
    max_triangles = 2 * n - 5

    # arrays that will store the triangulation graph
    triangles = List.duplicate(0, max_triangles * 3)
    halfedges = List.duplicate(0, max_triangles * 3)

    #temporary arrays for tracking the edges of the advancing convex hull
    hash_size = :math.sqrt(n)
                |> Float.ceil
                |> trunc
    hull_prev = List.duplicate(0, n) # edge to prev edge
    hull_next = List.duplicate(0, n) # edge to next edge
    hull_tri = List.duplicate(0, n) # edge to adjacent triangle
    hull_hash = List.duplicate(-1, hash_size) # angular edge hash

    {
      ids,
      _,
      min_x,
      min_y,
      max_x,
      max_y
    } =
      0..(n - 1)
      |> Enum.to_list
      |> List.foldl(
           {List.duplicate(0, n), 0, nil, nil, nil, nil},
           fn _, {
             ids,
             i,
             min_x,
             min_y,
             max_x,
             max_y
           } ->

             x = coords
                 |> Enum.at(2 * i)
             y = coords
                 |> Enum.at(2 * i + 1)
             min_x = if (min_x == nil || x < min_x), do: x, else: min_x
             min_y = if (min_y == nil || y < min_y), do: y, else: min_y
             max_x = if (max_x == nil || x > max_x), do: x, else: max_x
             max_y = if (max_y == nil || y > max_y), do: y, else: max_y
             {
               ids
               |> List.replace_at(i, i),
               i + 1,
               min_x,
               min_y,
               max_x,
               max_y
             }
           end
         )

    cx = (min_x + max_x) / 2
    cy = (min_y + max_y) / 2

    {_, _, i0} = findCenter(coords, n, cx, cy)
    i0x = coords
          |> Enum.at(2 * i0)
    i0y = coords
          |> Enum.at(2 * i0 + 1)

    {_, _, i1} = findClosestCenterPoint(coords, n, i0, i0x, i0y)
    i1x = coords
          |> Enum.at(2 * i1)
    i1y = coords
          |> Enum.at(2 * i1 + 1)

    {_, min_radius, i2} = findSmallestCircle(coords, n, i0, i0x, i0y, i1, i1x, i1y)
    i2x = coords
          |> Enum.at(2 * i2)
    i2y = coords
          |> Enum.at(2 * i2 + 1)

    if (min_radius == nil) do
      raise "No Delaunay triangulation exists for this input."
    end

    # swap the order of the seed points for counter-clockwise orientation
    {i1, i1x, i1y, i2, i2x, i2y} =
      if (Utils.orient(i0x, i0y, i1x, i1y, i2x, i2y)) do
        i = i1
        x = i1x
        y = i1y
        {i2, i2x, i2y, i, x, y}
      else
        {i1, i1x, i1y, i2, i2x, i2y}
      end

    {cx, cy} = Utils.circumcenter(i0x, i0y, i1x, i1y, i2x, i2y)

    dists =
      0..(n - 1)
      |> Enum.to_list
      |> List.foldl(
           List.duplicate(0, n),
           fn i, l ->
             l
             |> List.replace_at(
                  i,
                  Utils.dist(
                    coords
                    |> Enum.at(2 * i),
                    coords
                    |> Enum.at(2 * i + 1),
                    cx,
                    cy
                  )
                )
           end
         )

    # sort the points by distance from the seed triangle circumcenter
    ids = Utils.quicksort(ids, dists, 0, n - 1)

    # set up the seed triangle as the starting hull
    hull_start = i0
    hull_size = 3

    delaunay = %Delaunay{
      coords: coords,
      triangles: triangles,
      halfedges: halfedges,
      hash_size: hash_size,
      hull_start: hull_start,
      hull_prev: hull_prev,
      hull_next: hull_next,
      hull_tri: hull_tri,
      hull_hash: hull_hash,
      hull_size: hull_size,
      cx: cx,
      cy: cy,
      edge_stack: List.duplicate(0, 512)
    }

    hull_next = hull_next
                |> List.replace_at(i0, i1)
                |> List.replace_at(i1, i2)
                |> List.replace_at(i2, i0)

    hull_prev = hull_prev
                |> List.replace_at(i2, i1)
                |> List.replace_at(i0, i2)
                |> List.replace_at(i1, i0)

    hull_tri = hull_tri
               |> List.replace_at(i0, 0)
               |> List.replace_at(i1, 1)
               |> List.replace_at(i2, 2)

    hull_hash = hull_hash
                |> List.replace_at(Internals.hashKey(delaunay, i0x, i0y), i0)
                |> List.replace_at(Internals.hashKey(delaunay, i1x, i1y), i1)
                |> List.replace_at(Internals.hashKey(delaunay, i2x, i2y), i2)

    new_delaunay = %{
                     delaunay
                   |
                     triangles_len: 0,
                     hull_start: hull_start,
                     hull_prev: hull_prev,
                     hull_next: hull_next,
                     hull_tri: hull_tri,
                     hull_hash: hull_hash
                   }
                   |> Internals.addTriangle(i0, i1, i2, -1, -1, -1)
                   |> elem(0)

    new_delaunay =
      0..(length(ids) - 1)
      |> Enum.reduce(
           {
             new_delaunay,
             0,
             0,
             ids,
             coords,
             i0,
             i1,
             i2
           },
           &addTriangles/2
         )
      |> elem(0)
      |> (fn %Delaunay{hull_size: hull_size} = d -> %{d | hull: List.duplicate(0, hull_size)} end).()

    0..(new_delaunay.hull_size - 1)
    |> Enum.reduce(
         {new_delaunay, nil},
         &fillHull/2
       )
    |> elem(0)
    |> (
         fn %Delaunay{triangles: triangles, halfedges: halfedges, triangles_len: triangles_len} = d ->
           %{
             d
           |
             # get rid of temporary arrays
             hull_prev: nil,
             hull_next: nil,
             hull_tri: nil,
             # trim typed triangle mesh arrays
             triangles: triangles
                        |> Enum.slice(0, triangles_len),
             halfedges: halfedges
                        |> Enum.slice(0, triangles_len)
           }
         end).()
  end

  @doc """

  """
  def getTriangles(
        %Delaunay{
          triangles: triangles
        } = delaunay,
        i \\ 0
      ) do
    if (i < length(triangles)) do
      p0 = triangles
           |> Enum.at(i)
      p1 = triangles
           |> Enum.at(i + 1)
      p2 = triangles
           |> Enum.at(i + 2)

      [
        [
          getPoint(delaunay, p0),
          getPoint(delaunay, p1),
          getPoint(delaunay, p2)
        ]
        | getTriangles(delaunay, i + 3)
      ]
    else
      []
    end
  end

  @doc """

  """
  def getEdgePairs(
        %Delaunay{
          triangles: triangles
        } = delaunay,
        i \\ 0
      ) do
    if (i < length(triangles)) do
      p0 = triangles
           |> Enum.at(i)
      p1 = triangles
           |> Enum.at(i + 1)
      p2 = triangles
           |> Enum.at(i + 2)

      [
        (if (p0 < p1), do: {p0, p1}, else: {p1, p0}),
        (if (p1 < p2), do: {p1, p2}, else: {p2, p1}),
        (if (p0 < p2), do: {p0, p2}, else: {p2, p0})
        | getEdgePairs(delaunay, i + 3)
      ]
    else
      []
    end
  end

  @doc """

  """
  def getEdgePairsWithWeigth(delaunay) do
    getEdgePairs(delaunay)
    |> Enum.uniq
    |> Enum.map(
         fn {p0, p1} ->
           {p0x, p0y} = getPoint(delaunay, p0)
           {p1x, p1y} = getPoint(delaunay, p1)

           {p0, p1, Utils.dist(p0x, p0y, p1x, p1y)}
         end
       )
  end

  defp getPoint(%Delaunay{coords: coords}, p) do
    px = coords
         |> Enum.at(2 * p)
    py = coords
         |> Enum.at(2 * p + 1)
    {px, py}
  end


  defp addTriangles(
         k,
         {
           %Delaunay{
             hash_size: hash_size,
             hull_prev: hull_prev,
             hull_next: hull_next,
             hull_tri: hull_tri,
             hull_hash: hull_hash,
             hull_size: hull_size
           } = delaunay,
           xp,
           yp,
           ids,
           coords,
           i0,
           i1,
           i2
         } = acc
       ) do
    i = ids
        |> Enum.at(k)
    x = coords
        |> Enum.at(2 * i)
    y = coords
        |> Enum.at(2 * i + 1)

    cond do
      # skip near-duplicate points
      (k > 0 && abs(x - xp) <= @epsilon && abs(y - yp) <= @epsilon) ->
        acc
      # skip seed triangle points
      (i == i0 || i == i1 || i == i2) ->
        acc
      true ->
        # find a visible edge on the convex hull using edge hash
        key = delaunay
              |> Internals.hashKey(x, y)
        start =
          0..(hash_size - 1)
          |> Enum.reduce_while(
               0,
               fn j, _ ->
                 hull_k = rem(key + j, hash_size)
                 start = hull_hash
                         |> Enum.at(hull_k)
                 if (
                      start != -1 && start != hull_next
                                              |> Enum.at(start)) do
                   {:halt, start}
                 else
                   {:cont, start}
                 end
               end
             )
          |> (
               fn s ->
                 hull_prev
                 |> Enum.at(s)
               end).()

        e = findClosePoint(coords, hull_next, x, y, start, start)
        if (e == -1) do
          # likely a near-duplicate point; skip it
          acc
        else
          # add the first triangle from the point
          {new_delaunay, t} = Internals.addTriangle(
            delaunay,
            e,
            i,
            hull_next
            |> Enum.at(e),
            -1,
            -1,
            hull_tri
            |> Enum.at(e)
          )

          # recursively flip triangles from the point until they satisfy the Delaunay condition
          {new_delaunay, ar} = Internals.legalize(new_delaunay, t + 2)
          new_hull_tri = hull_tri
                         |> List.replace_at(i, ar)
                         |> List.replace_at(e, t) # keep track of boundary triangles on the hull

          new_delaunay = %{new_delaunay | hull_tri: new_hull_tri, hull_size: hull_size + 1}

          # walk forward through the hull, adding more triangles and flipping recursively
          n = hull_next
              |> Enum.at(e)

          {new_delaunay, n} = addTrianglesForward(
            new_delaunay,
            coords,
            x,
            y,
            n,
            i
          )

          # walk backward from the other side, adding more triangles and flipping
          {new_delaunay, e} =
            if (e != start) do
              {new_delaunay, e}
            else
              addTrianglesBackward(
                new_delaunay,
                coords,
                x,
                y,
                e,
                i
              )
            end

          # update the hull indices
          new_hull_prev = hull_prev
                          |> List.replace_at(i, e)
                          |> List.replace_at(n, i)
          new_hull_next = hull_next
                          |> List.replace_at(e, i)
                          |> List.replace_at(i, n)

          new_delaunay = %{
            new_delaunay
          |
            hull_start: e,
            hull_prev: new_hull_prev,
            hull_next: new_hull_next
          }

          # save the two new edges in the hash table
          new_hull_hash = hull_hash
                          |> List.replace_at(Internals.hashKey(new_delaunay, x, y), i)
                          |> List.replace_at(
                               Internals.hashKey(
                                 new_delaunay,
                                 coords
                                 |> Enum.at(2 * e),
                                 coords
                                 |> Enum.at(2 * e + 1)
                               ),
                               e
                             )

          {
            %{
              new_delaunay
            |
              hull_hash: new_hull_hash
            },
            xp,
            yp,
            ids,
            coords,
            i0,
            i1,
            i2
          }
        end
    end
  end

  defp addTrianglesForward(
         %Delaunay{hull_next: hull_next, hull_tri: hull_tri, hull_size: hull_size} = delaunay,
         coords,
         x,
         y,
         n,
         i
       ) do
    q = hull_next
        |> Enum.at(n)
    if (
         Utils.orient(
           x,
           y,
           coords
           |> Enum.at(2 * n),
           coords
           |> Enum.at(2 * n + 1),
           coords
           |> Enum.at(2 * q),
           coords
           |> Enum.at(2 * q + 1)
         )) do
      {new_delaunay, t} = Internals.addTriangle(
        delaunay,
        n,
        i,
        q,
        hull_tri
        |> Enum.at(i),
        -1,
        hull_tri
        |> Enum.at(n)
      )
      {new_delaunay, ar} = Internals.legalize(new_delaunay, t + 2)

      addTrianglesForward(
        %{
          new_delaunay
        |
          hull_next: hull_next
                     |> List.replace_at(n, n),
          hull_tri: hull_tri
                    |> List.replace_at(i, ar),
          hull_size: hull_size - 1,
        },
        coords,
        x,
        y,
        q,
        i
      )
    else
      {delaunay, n}
    end
  end

  defp addTrianglesBackward(
         %Delaunay{
           hull_next: hull_next,
           hull_prev: hull_prev,
           hull_tri: hull_tri,
           hull_size: hull_size
         } = delaunay,
         coords,
         x,
         y,
         e,
         i
       ) do
    q = hull_prev
        |> Enum.at(e)
    if (
         Utils.orient(
           x,
           y,
           coords
           |> Enum.at(2 * q),
           coords
           |> Enum.at(2 * q + 1),
           coords
           |> Enum.at(2 * e),
           coords
           |> Enum.at(2 * e + 1)
         )) do
      {new_delaunay, t} = Internals.addTriangle(
        delaunay,
        q,
        i,
        e,
        -1,
        hull_tri
        |> Enum.at(e),
        hull_tri
        |> Enum.at(q)
      )
      {new_delaunay, _ar} = Internals.legalize(new_delaunay, t + 2)

      addTrianglesBackward(
        %{
          new_delaunay
        |
          hull_next: hull_next
                     |> List.replace_at(e, e),
          hull_tri: hull_tri
                    |> List.replace_at(q, t),
          hull_size: hull_size - 1,
        },
        coords,
        x,
        y,
        q,
        i
      )
    else
      {delaunay, e}
    end
  end

  defp findClosePoint(coords, hull_next, x, y, start, e) do
    q = hull_next
        |> Enum.at(e)
    if (
         Utils.orient(
           x,
           y,
           coords
           |> Enum.at(2 * e),
           coords
           |> Enum.at(2 * e + 1),
           coords
           |> Enum.at(2 * q),
           coords
           |> Enum.at(2 * q + 1)
         )) do
      e
    else
      if (q == start) do
        -1
      else
        findClosePoint(coords, hull_next, x, y, start, q)
      end
    end
  end

  # pick a seed point close to the center
  defp findCenter(coords, n, cx, cy) do
    0..(n - 1)
    |> Enum.to_list
    |> List.foldl(
         {0, nil, nil},
         fn _, {i, min_dist, i0} ->
           d = Utils.dist(
             cx,
             cy,
             coords
             |> Enum.at(2 * i),
             coords
             |> Enum.at(2 * i + 1)
           )
           if (min_dist == nil || d < min_dist) do
             {i + 1, d, i}
           else
             {i + 1, min_dist, i0}
           end
         end
       )
  end

  # find the point closest to the seed
  defp findClosestCenterPoint(coords, n, i0, i0x, i0y) do
    0..(n - 1)
    |> Enum.to_list
    |> List.foldl(
         {0, nil, nil},
         fn _, {i, min_dist, i1} ->
           if (i == i0) do
             {i + 1, min_dist, i1}
           else
             d = Utils.dist(
               i0x,
               i0y,
               coords
               |> Enum.at(2 * i),
               coords
               |> Enum.at(2 * i + 1)
             )
             if ((min_dist == nil || d < min_dist && d > 0)) do
               {i + 1, d, i}
             else
               {i + 1, min_dist, i1}
             end
           end
         end
       )
  end

  # find the third point which forms the smallest circumcircle with the first two
  defp findSmallestCircle(coords, n, i0, i0x, i0y, i1, i1x, i1y) do
    0..(n - 1)
    |> Enum.to_list
    |> List.foldl(
         {0, nil, nil},
         fn _, {i, min_radius, i2} ->
           if (i == i0 || i == i1) do
             {i + 1, min_radius, i2}
           else
             r = Utils.circumradius(
               i0x,
               i0y,
               i1x,
               i1y,
               coords
               |> Enum.at(2 * i),
               coords
               |> Enum.at(2 * i + 1)
             )
             if (min_radius == nil || r < min_radius) do
               {i + 1, r, i}
             else
               {i + 1, min_radius, i2}
             end
           end
         end
       )
  end

  defp fillHull(
         i,
         {
           %Delaunay{
             hull: hull,
             hull_start: hull_start,
             hull_next: hull_next
           } = delaunay,
           e
         }
       ) do
    e = if (e != nil), do: e, else: hull_start
    {
      %{
        delaunay |
        hull: hull
              |> List.replace_at(i, e)
      },
      hull_next
      |> Enum.at(e)
    }
  end
end