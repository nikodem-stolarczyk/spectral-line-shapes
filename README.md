# spectral-line-shapes

The spectral-line-shapes repository hosts the implementations of the modified Hartman-Tran (`mHT`) spectral line-shape profile. The mHT is a variation of the Hartman-Tran profile (HT or HTP, introduced in 10.1016/j.jqsrt.2013.06.015). The article describing the mHT spectral lineshape function is currently under preparation. The mHT is a function with complex output. Its real and imaginary parts are the absorbtion and dispersion profiles, respectively.

The spectral-line-shapes repository features the implementations of the complex problability function (cpf), which are used to calculate the mHT. Though optimized specifically to work with the mHT, the cpf functions can also be used separately. The repository features two cpf functions: `cpf_accurate` (maintaining high accuracy, based on highly-modified 64-term Weideman algorithm jstor.org/stable/2158232), and `cpf_fast` (faster but less accurate, utilizing Humlicek's aproach in its first subregion, 10.1016/0022-4073(82)90078-4, and 24-term Weideman approach elsewhere jstor.org/stable/2158232). The detailed description of the modifications of the cpf functions will be summarized in the upcoming article.

At the moment, the python, fortran and matlab implementations have been developed, while the mathematica and labview implementations are under preparation.

## Installation

Clone the repository using:

```
git clone https://github.com/nikodem-stolarczyk/spectral-line-shapes/
```
This will clone all the implementations.

Alternatively, one can clone only one implementation. For instance, if one desires to clone the matlab implementation, they should use:

```
git clone https://github.com/nikodem-stolarczyk/spectral-line-shapes/matlab/
```
 
## Dependencies

Depending on the implementation, the following program versions are required

### python
- Python 3.X
- numpy and matplotlib (used solely in example_plots.py, the mHT and cpf functions do not need these libraries)
### Fortran
- **Fortran Standard**: The code is written in Fortran 90 (compatible with later standards). It's important to use a compiler that supports at least Fortran 90.
- **Recommended Compilers**: 
  - [GNU Fortran Compiler (gfortran)](https://gcc.gnu.org/wiki/GFortran)
  - [Intel Fortran Compiler (ifort)](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/fortran-compiler.html)
- **Build System**:
  - [Make](https://www.gnu.org/software/make/)
### matlab

### mathematica

### labview

## Usage

### python
The `example_absorption.py` and `example_dispersion.py` files present the basic syntax and usage of the mHT function to calculate absorption and dispersion of a spectral line. The `example_mHT_optional_parameters.py` presents the syntax with the optional parameters.
The `example_plots.py` file generates two plots:
- the absorption and dispersion profiles of the Ar-perturbed 3-0 S(1) line in H<sub>2</sub> (see [Stolarczyk et al.](https://doi.org/10.1063/5.0139229))
![plot1](media/plot_example_1.svg)
- the comparison the absorption profiles of the Ar- and He-perturbed 3-0 S(1) line H<sub>2</sub> (see [Stolarczyk et al.](https://doi.org/10.1063/5.0139229) and [Slowinski et al.](https://doi.org/10.1103/PhysRevA.101.052705)).
![plot2](media/plot_example_2.svg)

The `example_cpf` file presents the syntax for the `cpf_accurate` and `cpf_fast` functions.

### Fortran

### matlab

### mathematica

### labview
## Data Format

### mHT 

The `mHT` function requires the following arguments:
- `nu0` - Unperturbed line position in cm-1.
- `GamD` - Doppler broadening in cm-1
- `Gam0` - Speed-averaged line-width in cm-1.
- `Gam2` - Speed dependence of the line-width in cm-1.
- `Shift0` - Speed-averaged line-shift in cm-1.
- `Shift2` - Speed dependence of the line-shift in cm-1.
- `NuOptRe` - Real part of the Dicke parameter in cm-1.
- `NuOptIm` - Imaginary part of the Dicke parameter in cm-1.
- `nu` - Current WaveNumber of the Computation in cm-1.

The following arguments are optional:
- `Ylm` - Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
- `Xlm` - Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
- `alpha` - Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless.

### cpf

Both `cpf_accurate` and `cpf_fast` functions require two arguments, `x` and `y`. The meaning of `x` and `y` input values varies, depending on the application of the cpf. In the spectral line-shape context, `y` is the ratio between Lorentian and Doppler broadening (y=sqrt(log(2))*Gam0/GamD), and `x` is the normalized detuning from the line center (x=sqrt(log(2))*(nu-nu0)/GamD). The `x` and `y` values are dimensionless.    

## Contributing

Contributions are welcome. Please fork the repository and submit pull requests with your enhancements.
