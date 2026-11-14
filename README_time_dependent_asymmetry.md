# Azimuthal Asymmetry and Rigid Rotation

The `oimComponentRadialProfile` implementation supports azimuthally modulated
radial brightness profiles, consistently in both image and Fourier space:

$$
I(r,\phi,\lambda,t)=I_0(r,\lambda,t)\left[
1+\sum_{m=1}^{M}a_m(\lambda,t)
\cos\left(m\left[\phi-\phi_m(\lambda,t)\right]\right)\right].
$$

## Asymmetry handling

- **Image consistency:** Azimuthal modulation is evaluated in the component's
  intrinsic (deprojected) coordinates in `getImage()`.
- **Fourier calculation:** Harmonic `m` contributes a `J_m` Bessel transform
  with the complex factor `(-1j)**m` and the corresponding azimuthal phase.
- **Observation-dependent parameters:** `skw1`, `skwPa1`, and higher-order
  amplitudes/phases are evaluated at each observation's wavelength and epoch,
  rather than being baked into the shared radial Fourier kernels. They can
  therefore use the existing wavelength/time parameter interpolators.
- **Multiple harmonics:** A harmonic order can be selected with `modulation`.
  Higher-order amplitudes default to zero, and parameter names are unique.
- **Backward compatibility:** With `asymmetric=False`, the symmetric
  calculation follows the existing implementation.

## Optional rigid rotation

Two optional parameters define a rigidly rotating azimuthal pattern:

- `rotPeriod`: signed rotation period in days. A negative value reverses
  the direction of rotation; zero is not allowed.
- `rotT0`: reference epoch in MJD. Both rotation parameters must be supplied
  together, and rotation requires `asymmetric=True`.

The phase of every harmonic evolves as

$$
\phi_m(t)=\phi_{m,0}+\frac{2\pi(t-t_0)}{P}.
$$

All harmonics share this same angular speed. This models **rigid pattern
rotation**, not Keplerian differential rotation. The component-level methods
use `rotT0` when no epoch is provided.

## Automated tests

The test module is `tests/test_radial_asymmetry.py`. Run it with:

```bash
python -m pytest tests/test_radial_asymmetry.py -v
```

The module contains 13 test functions (15 pytest cases, since two tests run
with two parameter values each). It checks:

1. Harmonic parameter initialization and invalid rotation configurations.
2. Image periodicity and positive/negative rigid-rotation directions.
3. Complex Fourier fluxes against independent direct polar integration for
   first-, second-, and third-order azimuthal harmonics.
4. Observation-by-observation wavelength/time interpolation, including
   intermediate values, using separately configured static reference models.
5. Hermitian Fourier symmetry and conservation of zero-frequency flux.
6. Symmetric-mode backward compatibility and an independent high-resolution
   radial Hankel integral, including the absolute flux normalization.
7. Consistency of `getImage()` pixels with `getComplexCoherentFlux()` through
   a direct discrete Fourier sum, with ellipticity, position angle, and offsets.
8. Parameter updates after initialization in both image and Fourier space.
9. Asymmetry and rotation of `oimTempGrad` in the L and N bands, comparing
   rotating and fixed-phase models in both image and Fourier space.

**Interpretation and limitations:** The original polar-integral tests use
independently constructed azimuthal brightness distributions but share the
model's radial quadrature grid/weights, and therefore compare *normalized*
complex fluxes. The independent symmetric radial-integral test also checks
the absolute flux scale. The pixel-based Fourier comparison is affected by
finite pixel size and field of view. Parameter-interpolation, parameter-update,
and `oimTempGrad` tests check that the model behaves as expected; they are not
independent checks of the Fourier formula.
The full existing package test suite should also be run before merging.
