Bidirectional Path Tracer
====================================================================================================

![BDPT Output 25 Samples][current25]

25 Samples

![BDPT Output 1000 Samples][current1000]

1000 Samples

Background
------------------
This BDPT implementation written by Isaac Kousari is built on top of the _Ray Tracing in One Weekend_ series of books.

### v1.0 Implementation details
  - Cosine Hemisphere sampling for lambertian reflections
  - Dieletrics with schlick approximation for refraction
  - Russian Roulette termination in camera and light path creation
  - Balanced Heuristic MIS

### TODOs

  - Fix fireflies
  - Fix dieletric reflections (_I think partially to do with reflection probabilty being assigned to 1_)
    - Light is not being reflected on top of sphere (_I think is a corner case in the delta material handling at vertex connection_)
  - Improve caustics
  - Make README.md nicer :)
  - Make easier to dump the individual unweighted s & t contributions
  - Finish looking into Microfacets and implement Disney Principled BSDF
  - Improve build system to proper .h/.cpp files

Directory Structure
-------------------
  - `images/` --
    Contains all of the current results

  - `src/` --
    Contains the source.
    BDPT implementation is inside `bdptcamera.h` because the original authors wrote a very very hacky build script v1.1 or v1.2 should move this to proper class organization and fix this.

Building and Running
---------------------
Copies of the source are provided for you to check your work and compare against. If you wish to
build the provided source, this project uses CMake. To build, go to the root of the project
directory and run the following commands to create the debug version of every executable:

    $ cmake -B build
    $ cmake --build build

You should run `cmake -B build` whenever you change your project `CMakeLists.txt` file (like when
adding a new source file).

### Optimized Builds
CMake supports Release and Debug configurations. These require slightly different invocations
across Windows (MSVC) and Linux/macOS (using GCC or Clang). The following instructions will place
optimized binaries under `build/Release` and debug binaries (unoptimized and containing debug
symbols) under `build/Debug`:

On Windows:

```shell
$ cmake -B build
$ cmake --build build --config Release  # Create release binaries in `build\Release`
$ cmake --build build --config Debug    # Create debug binaries in `build\Debug`
```

On Linux / macOS:

```shell
# Configure and build release binaries under `build/Release`
$ cmake -B build/Release -DCMAKE_BUILD_TYPE=Release
$ cmake --build build/Release

# Configure and build debug binaries under `build/Debug`
$ cmake -B build/Debug -DCMAKE_BUILD_TYPE=Debug
$ cmake --build build/Debug
```

[current25]:          images/bdpt.png
[current1000]:          images/bdpt1000.png
