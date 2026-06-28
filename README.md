# RY-SlopeSearch
## Introduction
Here is the source code of RY-SlopeSearch, which is an slope-based alignment-free algorithms for fast DNA seqs comparision.

All ranking results of RY-SlopeSearch are publicly available on the AFproject website. You can find them by searching for the keyword `RY-SlopeSearch-v1.0.0` on the AFproject website.

### $F(k)$ Curve for Word Lengths $k$

We used two sequences, *Shigella dysenteriae* Sd197 and *Escherichia coli* UTI89, to generate the $F(k)$ curve for word lengths $k$. This curve provides a visualization of $F(k)$ and helps illustrate the slope-based method.

![Relationship between F(k) and word length k](./doc/f_k.png)

**Figure:** Relationship between $F(k)$ and word length $k$ for $k \in [2, 24]$, computed for *Shigella dysenteriae* Sd197 and *Escherichia coli* UTI89. The curve was generated using all $k$-mers and many-to-many matching.

### Empirical Determination of the $k$-Range for Slope Estimation

The slope of $F(k)$ is estimated over a restricted interval of $k$ values. Selecting an appropriate range is important: for small $k$, background matches dominate and introduce noise, whereas for large $k$, the number of observed matches becomes sparse.

In our implementation, the interval is determined by an empirical formula based on the average sequence length, following Slope-SpaM.

Let

$$
L_{\mathrm{avg}} = \frac{L_1 + L_2}{2}
$$

denote the average length of the two sequences, and let $\ell$ denote the pattern length. For full $k$-mers, $\ell = 1$.

We define

$$
k_{\min}=\max\left(\left\lceil\frac{\log(L_{\mathrm{avg}}) + 0.69}{0.875}\right\rceil,\ell\right)
$$

and

$$
k_{\max}=\max\left(\left\lfloor\frac{\log(L_{\mathrm{avg}})}{0.634}\right\rfloor,\ell\right).
$$

These constants are taken from the Slope-SpaM paper. The slope is then fitted using all integer values

$$
k \in [k_{\min}, k_{\max}].
$$

## Getting Start
### Requirements
- GCC 12.3.0 or higher
- CMake: Version 3.20 or higher
- third party
  - SeqAn3 Library
  - range-v3 Library: Needed as a dependency for SeqAn3
  - yaml-cpp Library: Used for YAML configuration parsing
  - OpenMP: Optional, but if enabled, it requires a compatible OpenMP installation for parallel processing support

### Installation
After downloading the code, run following code:
```shell
cd RY-SlopeSearch
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc-12 -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 ..
make
```

*ignore `warning: use of ‘std::hardware_destructive_interference_size’`*


### How to run an example
#### Example location
- Example Folder: `RY-SlopeSearch/example`
- Example Data: `RY-SlopeSearch/example/example_data/assembled-ecoli`

#### Step 1: Create an output folder
- Example Output Folder: `RY-SlopeSearch/example/example_output`

#### Step 2: Editing configure
config location: `RY-SlopeSearch/example/example_config.yaml`
```yaml
dataset:
  datasets_directory: "(RY-SlopeSearch Project Root)/example/example_data/assembled-ecoli"

parameters:
  using_one_to_one_matching: true
  using_background_matches: false
  use_openmp: true

# methods
#   basic_kmer_matches
#   start_ry_matches
#   start_ry_4_6_matches
#   start_ry_4_9_matches
#   start_ry_16_matches
#   start_ry_32_matches
#   start_ry_64_matches
#   start_ry_128_matches
sampling_method: "start_ry_128_matches"

output_directory: "(RY-SlopeSearch Project Root)/example/example_output"

use_progress_bar: true
```
The parameters of configs

| Parameter                     | Description                                                                          |
|-------------------------------|--------------------------------------------------------------------------------------|
| `dataset.datasets_directory`  | Specifies the directory path for the sequence data to be used.                       |
| `parameters.using_one_to_one_matching` | Enables the one-to-one matching algorithm if set to true.                |
| `parameters.using_background_matches`  | Enables background matches if set to true.                               |
| `parameters.use_openmp`       | Enables multi-threaded computation with OpenMP if supported in your environment.     |
| `sampling_method`             | Specifies the name of the sampling algorithm to use.                                  |
| `output_directory`            | Specifies the path where the output files will be saved.                             |
| `use_progress_bar`            | Enables progress bar when doing calculation if set to true.                       |

- `fasta` files should be collected under `dataset.datasets_directory`
- each `fasta` file has one genome
- All supported word sets are listed in the comments.
  - basic_kmer_matches: all k-mers
  - start_ry_matches: RY word set
  - start_ry_4_6_matches: RY 4-6 word set
  - start_ry_4_9_matches: RY 4-9 word set
  - start_ry_16_matches: RY16-11 word set
  - start_ry_32_matches: RY32-11 word set
  - start_ry_64_matches: RY64-13 word set
  - start_ry_128_matches: RY128-12 word set
- output:
  - a time stamp folder will be created under `output_directory`, your results are in there
  - results includes:
    - config you used
    - runtime log
    - `.phy` distance matrix
    - performance log
      - time of matrix computation

*All paths are required to be **absolute paths** .*

#### Step 3: Run Program
- Enter into build folder: `cd ./build`
- Run the program: `./RY-SlopeSearch --config (RY-SlopeSearch project root)/configs/example_config.yaml`
  - `--config`: the location of your config

#### Step 4: Check the result
All results are in `output_directory`

#### Step 5: Upload `.phy` file to AFProject website
You can upload your `.phy` file of fish_mito genomes (example data) to [AFProject-benchmark-fish_mito](https://afproject.org/app/benchmark/genome/std/assembled/fish_mito/) to check nRF/nQD and your ranking.

## Contact
If you have any questiones, please feel free to leave a message to me!
