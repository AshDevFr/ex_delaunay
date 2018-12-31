# Delaunay

Elixir library for Delaunay triangulation of 2D points

## Example

```elixir
points = [{111, 0}, {98, 77}, {57, 194}, {120, 110}, {28, 187}, {113, 61}, {43, 36}, ...]

delaunay =
  points  
  |> Delaunay.from

IO.inspect delaunay.triangles
// [1, 18, 5, 1, 3, 18, 1, 29, 3, 3, 29, 18, 15, 31, 18, 18, 31, 5,
   31, 14, 5, 31, 16, 14, 25, 17, 29, 29, 10, 18, 1, 25, 29, 8, 28, 5, ...],
```

## Installation

If [available in Hex](https://hex.pm/docs/publish), the package can be installed
by adding `delaunay` to your list of dependencies in `mix.exs`:

```elixir
def deps do
  [
    {:delaunay, "~> 0.1.0"}
  ]
end
```

Documentation can be generated with [ExDoc](https://github.com/elixir-lang/ex_doc)
and published on [HexDocs](https://hexdocs.pm). Once published, the docs can
be found at [https://hexdocs.pm/delaunay](https://hexdocs.pm/delaunay).
