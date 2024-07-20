# spectral-line-shapes

The spectral-line-shapes repository hosts the implementations of the modified Hartman-Tran (`mHT`) spectral line-shape profile. The mHT is a variation of the Hartman-Tran profile (HT or HTP, introduced in 10.1016/j.jqsrt.2013.06.015). The article describing the mHT spectral lineshape function is currently under preparation.

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
### fortran
- [Make](https://www.gnu.org/software/make/)
- A modern Fortran compiler, such as [GNU Fortran Compiler](https://gcc.gnu.org/wiki/GFortran) (gfortran) or
 [Intel Fortran Compiler](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/fortran-compiler.html) (ifort)
### matlab

### mathematica

### labview

## Usage

The `usage_example.x` (x stands for the extension depending on the implementation) provides examples of using the mHT, cpf_accurate and cpf_fast functions in each of the programing languages.

## Data Format

### mHT 

The `mHT` function requires the following arguments:
- `nu0` - Unperturbed line position in cm-1.
- `Gam0` - Speed-averaged line-width in cm-1.
- `Gam2` - Speed dependence of the line-width in cm-1.
- `Shift0` - Speed-averaged line-shift in cm-1.
- `Shift2` - Speed dependence of the line-shift in cm-1.
- `NuOptRe` - Real part of the Dicke parameter in cm-1.
- `NuOptIm` - Imaginary part of the Dicke parameter in cm-1.
- `nu` - Current WaveNumber of the Computation in cm-1.

The following arguments are optional:
- `Sw` - Statistical weight.
- `Ylm` - Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
- `Xlm` - Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
- `alpha` - Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless.

### cpf

Both `cpf_accurate` and `cpf_fast` functions require two arguments, `x` and `y`. The meaning of `x` and `y` input values varies, depending on the application of the cpf. In the spectral line-shape context, `y` is the ratio between Lorentian and Doppler broadening (y=sqrt(log(2))*Gamma_L/Gamma_D), and `x` is the normalized detuning from the line center (x=sqrt(log(2))*(nu-nu_0)/Gamma_D). The `x` and `y` values are dimensionless.    

## Contributing

Contributions are welcome. Please fork the repository and submit pull requests with your enhancements.
