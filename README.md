# Cornea detection

Julia software for detecting the cornea in a fluorescence image. This was developed and run in Julia 1.5.3, with the packages and versions listed in the `Manifest.toml` file.

## Installation

See the [documentation on the package manager](https://julialang.github.io/Pkg.jl/v1/) if you need additional information.

1. (optional) Create a directory, which I will designate here as `foo`.
2. Start Julia.
3. Enter `using Pkg`.
4. (recommended) Enter `pkg"activate foo"`. This creates a separate environment for the software, to help avoid dependency version conflicts.
5. Enter `pkg"add https://github.com/tobydriscoll/cornea-detection"`. This will download all required dependencies, which will take a few minutes.

## Usage

If you created a dedicated environment for installation, you will have to follow the process above to activate it in each new Julia session. Then enter

```julia
using CorneaDetection
```

On the first invocation it will take several minutes to compile all the dependencies. Loading should be much faster after that.

The fundamental entry point is the `detect` method. Running `detect(img)` for image `img` performs the cornea detection and returns a named tuple with these primary fields:

* `cenrow` row location of the cornea center, normalized by the image height
* `cencol` column location of the cornea center, normalized by the image height
* `radius` radius of the cornea center, normalized by the image height

There are additional diagnostic fields you can ignore.

If you run `detect(folder)` with `folder` being a vector of file names, the images will be processed serially. In this mode, the final state of one image is used as a candidate initial state for the next one, making a continuous transition more likely for the optimizer.

The file `src/parameters.jl` defines parameter values that affect the behavior of the algorithm. Some may need to be adjusted to get the most consistent results on a new data set. The most important ones may be those related to automatic detection of the Purkinje image, which is used to find candidates to initialize the optimization.

## Sample data

You can find [sample data](https://doi.org/10.5281/zenodo.4426013) from two videos to try out the algorithms.
