# geometry

## Build

Requires CMake 3.20+ and a C++17 compiler.

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
