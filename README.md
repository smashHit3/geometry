# Geometry Visualizer

An interactive Qt Widgets application for creating and manipulating simple geometry.

## Build

Requires CMake 3.20+, a C++17 compiler, and Qt 5 or Qt 6 with the Widgets component.

```sh
cmake -S . -B build
cmake --build build
```

## Test

```sh
ctest --test-dir build --output-on-failure
```

Unit tests for individual modules are located in each module's `unittest/` directory.

## Run

```sh
./build/geometry
```

Use the toolbar or shortcuts to add rectangles and circles. Select a shape and drag it
to move it; use `Delete` to remove selected shapes.
