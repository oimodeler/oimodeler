"""Tests for radial-profile asymmetry and rotation."""

import numpy as np
import pytest
from scipy.integrate import simpson
from scipy.special import j0
import oimodeler as oim


class GaussianTestRing(oim.oimComponentRadialProfile):
    """Gaussian ring used for the tests."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._r = np.linspace(0.01, 5.0, 256)

    def _radialProfileFunction(self, r, wl, t):
        return np.exp(-0.5 * ((r - 2.0) / 0.4)**2)


def test_asymmetry_initialization():
    """Check skew parameters and default values."""

    ring = GaussianTestRing(
        asymmetric=True,
        modulation=2,
        skw1=0.3,
        skwPa1=30.0,
    )

    assert ring.asymmetric
    assert ring.modulation == 2

    assert ring.skw1.value == 0.3
    assert ring.skwPa1.value == 30.0

    assert ring.skw2.value == 0.0

    assert ring.skw1.name == "skw1"
    assert ring.skwPa1.name == "skwPa1"


def test_rigid_rotation():
    """Check that the image repeats after one rotation."""

    ring = GaussianTestRing(
        asymmetric=True,
        modulation=1,
        skw1=0.4,
        skwPa1=30.0,
        rotPeriod=20.0,
        rotT0=60000.0,
    )

    images = ring.getImage(
        dim=80,
        pixSize=0.1,
        wl=np.array([3.5e-6]),
        t=np.array([60000.0, 60010.0, 60020.0]),
    )

    # Same image after one period
    np.testing.assert_allclose(
        images[0],
        images[2],
        rtol=0,
        atol=1e-12,
    )

    # Different image after half a period
    assert not np.allclose(images[0], images[1])


def test_invalid_rotation_parameters():
    """Check that invalid rotation settings give an error."""

    with pytest.raises(ValueError):
        GaussianTestRing(
            asymmetric=True,
            rotPeriod=0.0,
            rotT0=60000.0,
        )

    with pytest.raises(ValueError):
        GaussianTestRing(
            asymmetric=True,
            rotPeriod=20.0,
        )

    with pytest.raises(ValueError):
        GaussianTestRing(
            asymmetric=False,
            rotPeriod=20.0,
            rotT0=60000.0,
        )

def test_asymmetric_fourier_against_direct_integration():
    """Compare the first harmonic with a direct Fourier integral."""

    ring = GaussianTestRing(
        asymmetric=True,
        modulation=1,
        skw1=0.4,
        skwPa1=30.0,
        rotPeriod=20.0,
        rotT0=60000.0,
    )

    # Spatial frequencies in cycles/rad
    ucoord = np.array([0.0, 1e7, 2e7, -3e7, 4e7])
    vcoord = np.array([0.0, 2e7, -1e7, 2e7, -3e7])

    wavelength = 3.5e-6
    epoch = 60005.0

    wl = np.full(len(ucoord), wavelength)
    t = np.full(len(ucoord), epoch)

    # OImodeler result
    native = ring.getComplexCoherentFlux(ucoord, vcoord, wl, t)

    # Direct Fourier integral
    mas2rad = np.pi / (180.0 * 3600.0 * 1000.0)

    r_mas = ring.r
    r_rad = r_mas * mas2rad
    dr_rad = ring.dr * mas2rad

    phi = np.linspace(0, 2*np.pi, 720, endpoint=False)
    dphi = 2*np.pi / len(phi)

    # Gaussian profile
    radial = np.exp(-0.5 * ((r_mas - 2.0) / 0.4)**2)

    # Skew term
    phase = np.deg2rad(30.0) + 2*np.pi*(epoch - 60000.0)/20.0
    modulation = 1.0 + 0.4*np.cos(phi - phase)

    reference = []

    for uu, vv in zip(ucoord, vcoord):

        exponent = -2j*np.pi*r_rad[:, None] * (
            uu*np.cos(phi)[None, :] +
            vv*np.sin(phi)[None, :]
        )

        integrand = (
            radial[:, None]
            * modulation[None, :]
            * np.exp(exponent)
            * r_rad[:, None]
            * dr_rad[:, None]
        )

        reference.append(np.sum(integrand) * dphi)

    reference = np.asarray(reference)

    # Compare normalized results
    native = native / native[0]
    reference = reference / reference[0]

    np.testing.assert_allclose(
        native,
        reference,
        rtol=1e-4,
        atol=1e-6,
    )

@pytest.mark.parametrize("order", [2, 3])
def test_multiple_harmonics_fourier(order):
    """Compare higher harmonics with a direct Fourier integral."""

    amplitudes = [0.25, 0.12, 0.08]
    phases = [30.0, 75.0, -40.0]

    params = {
        "asymmetric": True,
        "modulation": order,
    }

    for m in range(1, order + 1):
        params[f"skw{m}"] = amplitudes[m - 1]
        params[f"skwPa{m}"] = phases[m - 1]

    ring = GaussianTestRing(**params)

    ucoord = np.array([0.0, 1e7, 2e7, -3e7, 4e7])
    vcoord = np.array([0.0, 2e7, -1e7, 2e7, -3e7])

    wl = np.full(len(ucoord), 3.5e-6)
    t = np.full(len(ucoord), 60000.0)

    native = ring.getComplexCoherentFlux(ucoord, vcoord, wl, t)

    # Direct integral
    mas2rad = np.pi / (180.0 * 3600.0 * 1000.0)

    r_mas = ring.r
    r_rad = r_mas * mas2rad
    dr_rad = ring.dr * mas2rad

    phi = np.linspace(0, 2*np.pi, 1024, endpoint=False)
    dphi = 2*np.pi / len(phi)

    radial = np.exp(-0.5 * ((r_mas - 2.0) / 0.4)**2)

    modulation = np.ones_like(phi)

    for m in range(1, order + 1):
        modulation += amplitudes[m-1] * np.cos(
            m * (phi - np.deg2rad(phases[m-1]))
        )

    reference = []

    for uu, vv in zip(ucoord, vcoord):

        exponent = -2j*np.pi*r_rad[:, None] * (
            uu*np.cos(phi)[None, :]
            + vv*np.sin(phi)[None, :]
        )

        integrand = (
            radial[:, None]
            * modulation[None, :]
            * np.exp(exponent)
            * r_rad[:, None]
            * dr_rad[:, None]
        )

        reference.append(np.sum(integrand) * dphi)

    reference = np.asarray(reference)

    native /= native[0]
    reference /= reference[0]

    np.testing.assert_allclose(
        native,
        reference,
        rtol=1e-4,
        atol=1e-6,
    )

def test_wavelength_and_time_dependent_asymmetry():
    """Check skew parameters that change with wavelength and time."""

    wavelengths = np.array([3.0e-6, 4.0e-6])
    epochs = np.array([60000.0, 60010.0])

    amplitudes = [0.1, 0.4]
    phases = [20.0, 100.0]

    ring = GaussianTestRing(
        asymmetric=True,
        modulation=1,
        skw1=oim.oimInterp(
            "wl", wl=wavelengths, values=amplitudes
        ),
        skwPa1=oim.oimInterp(
            "time", mjd=epochs, values=phases
        ),
    )

    ucoord = np.array([1e7, 2e7, -3e7, 4e7, -2e7, 3e7])
    vcoord = np.array([2e7, -1e7, 2e7, -3e7, 4e7, -4e7])

    # Include points between the given values
    wl = np.array([3e-6, 4e-6, 3e-6, 4e-6, 3.5e-6, 3.5e-6])
    t = np.array([60000., 60000., 60010., 60010., 60005., 60000.])

    # All points at once
    combined = ring.getComplexCoherentFlux(
        ucoord, vcoord, wl, t
    )

    # One model for each point
    reference = []

    for uu, vv, wli, ti in zip(ucoord, vcoord, wl, t):

        amp = np.interp(wli, wavelengths, amplitudes)
        phase = np.interp(ti, epochs, phases)

        static_ring = GaussianTestRing(
            asymmetric=True,
            modulation=1,
            skw1=float(amp),
            skwPa1=float(phase),
        )

        flux = static_ring.getComplexCoherentFlux(
            np.array([uu]),
            np.array([vv]),
            np.array([wli]),
            np.array([ti]),
        )

        reference.append(flux[0])

    np.testing.assert_allclose(
        combined,
        np.asarray(reference),
        rtol=1e-10,
        atol=1e-10,
    )

def test_hermitian_symmetry():
    """Check the Fourier values at opposite spatial frequencies."""

    ring = GaussianTestRing(
        asymmetric=True,
        modulation=2,
        skw1=0.3,
        skwPa1=45.0,
        skw2=0.2,
        skwPa2=-30.0,
        rotPeriod=20.0,
        rotT0=60000.0,
    )

    u = np.array([1e7, 2e7, -3e7])
    v = np.array([2e7, -1e7, 4e7])
    wl = np.full(3, 3.5e-6)
    t = np.full(3, 60005.0)

    positive = ring.getComplexCoherentFlux(u, v, wl, t)
    negative = ring.getComplexCoherentFlux(-u, -v, wl, t)

    np.testing.assert_allclose(
        negative,
        np.conj(positive),
        rtol=1e-12,
        atol=1e-12,
    )


def test_zero_baseline_flux():
    """Check that skew does not change the total flux."""

    symmetric = GaussianTestRing(asymmetric=False)

    asymmetric = GaussianTestRing(
        asymmetric=True,
        modulation=2,
        skw1=0.4,
        skwPa1=30.0,
        skw2=0.2,
        skwPa2=60.0,
    )

    u = np.array([0.0])
    v = np.array([0.0])
    wl = np.array([3.5e-6])
    t = np.array([60000.0])

    f_sym = symmetric.getComplexCoherentFlux(u, v, wl, t)
    f_asym = asymmetric.getComplexCoherentFlux(u, v, wl, t)

    np.testing.assert_allclose(
        f_sym, f_asym,
        rtol=1e-12,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    "period, expected_phase",
    [(20.0, 120.0), (-20.0, -60.0)],
)
def test_rotation_direction(period, expected_phase):
    """Check the direction of rotation for both period signs."""

    rotating = GaussianTestRing(
        asymmetric=True,
        skw1=0.4,
        skwPa1=30.0,
        rotPeriod=period,
        rotT0=60000.0,
    )

    static = GaussianTestRing(
        asymmetric=True,
        skw1=0.4,
        skwPa1=expected_phase,
    )

    u = np.array([1e7, 2e7, -3e7])
    v = np.array([2e7, -1e7, 4e7])
    wl = np.full(3, 3.5e-6)
    t = np.full(3, 60005.0)

    f_rot = rotating.getComplexCoherentFlux(u, v, wl, t)
    f_static = static.getComplexCoherentFlux(u, v, wl, t)

    np.testing.assert_allclose(
        f_rot, f_static,
        rtol=1e-12,
        atol=1e-12,
    )


def test_symmetric_fourier_regression():
    """Compare symmetric models with a direct radial integral."""
    symmetric = GaussianTestRing(asymmetric=False)
    zero_skew = GaussianTestRing(
        asymmetric=True, modulation=1, skw1=0.0, skwPa1=73.0,
    )

    u = np.array([0.0, 1e7, 2e7, -3e7, 4e7])
    v = np.array([0.0, 2e7, -1e7, 4e7, -3e7])
    wl = np.full(u.size, 3.5e-6)
    t = np.full(u.size, 60000.0)

    native = symmetric.getComplexCoherentFlux(u, v, wl, t)
    zero_amp = zero_skew.getComplexCoherentFlux(u, v, wl, t)
    np.testing.assert_allclose(zero_amp, native, rtol=1e-12, atol=1e-12)

    # Finer radial grid
    mas_to_rad = np.pi / (180.0 * 3600.0 * 1000.0)
    r_mas = np.linspace(0.0, 6.0, 6001)
    r_rad = r_mas * mas_to_rad
    intensity = np.exp(-0.5 * ((r_mas - 2.0) / 0.4)**2)
    spatial_frequency = np.hypot(u, v)
    integrand = (
        2.0 * np.pi * r_rad[:, None]
        * intensity[:, None]
        * j0(2.0 * np.pi * r_rad[:, None] * spatial_frequency[None, :])
    )
    # Convert to Jy
    reference = 1e23 * simpson(integrand, x=r_rad, axis=0)

    np.testing.assert_allclose(native.real, reference, rtol=2e-4, atol=1e-6)
    np.testing.assert_allclose(native.imag, 0.0, rtol=0, atol=1e-6)


def test_image_fourier_consistency():
    """Compare the Fourier result with the model image."""
    ring = GaussianTestRing(
        asymmetric=True,
        modulation=2,
        skw1=0.3,
        skwPa1=35.0,
        skw2=0.12,
        skwPa2=-28.0,
        rotPeriod=20.0,
        rotT0=60000.0,
        elong=1.25,
        pa=32.0,
        x=0.24,
        y=-0.18,
    )
    dim = 256
    pix_mas = 0.05
    wavelength = 3.5e-6
    epoch = 60005.0

    frame = ring.getImage(
        dim=dim,
        pixSize=pix_mas,
        wl=np.array([wavelength]),
        t=np.array([epoch]),
    )[0, 0]

    mas_to_rad = np.pi / (180.0 * 3600.0 * 1000.0)
    coords = (np.arange(dim) - dim // 2) * pix_mas * mas_to_rad
    xx, yy = np.meshgrid(coords, coords)
    x_flat, y_flat = xx.ravel(), yy.ravel()
    flux_flat = frame.ravel()

    u = np.array([0.0, 1e7, 2e7, -3e7, 4e7])
    v = np.array([0.0, 2e7, -1e7, 4e7, -3e7])
    wl = np.full(u.size, wavelength)
    t = np.full(u.size, epoch)

    direct_image_ft = np.array([
        np.sum(flux_flat * np.exp(-2j * np.pi * (uu*x_flat + vv*y_flat)))
        for uu, vv in zip(u, v)
    ])
    direct_image_ft /= np.sum(flux_flat)

    native = ring.getComplexCoherentFlux(u, v, wl, t)
    native /= native[0]

    np.testing.assert_allclose(
        direct_image_ft, native, rtol=2e-3, atol=5e-4,
    )


def test_parameter_updates_after_initialization():
    """Check the results after changing model parameters."""
    settings = dict(
        asymmetric=True, modulation=1,
        skw1=0.2, skwPa1=25.0,
        rotPeriod=20.0, rotT0=60000.0,
    )
    ring = GaussianTestRing(**settings)

    u = np.array([1e7, 2e7, -3e7])
    v = np.array([2e7, -1e7, 4e7])
    wl = np.full(u.size, 3.5e-6)
    t = np.full(u.size, 60005.0)

    def evaluate(component):
        fourier = component.getComplexCoherentFlux(u, v, wl, t)
        image = component.getImage(
            dim=96, pixSize=0.1,
            wl=np.array([3.5e-6]), t=np.array([60005.0]),
        )[0, 0]
        return fourier, image

    previous_flux, previous_image = evaluate(ring)

    for param, value in [("skw1", 0.5), ("skwPa1", 75.0), ("rotPeriod", -20.0)]:
        ring.params[param].value = value
        settings[param] = value
        fresh_ring = GaussianTestRing(**settings)

        flux, image = evaluate(ring)
        fresh_flux, fresh_image = evaluate(fresh_ring)

        np.testing.assert_allclose(flux, fresh_flux, rtol=1e-12, atol=1e-8)
        np.testing.assert_allclose(image, fresh_image, rtol=1e-12, atol=1e-12)
        assert not np.allclose(flux, previous_flux, rtol=1e-5, atol=1e-8)
        assert not np.allclose(image, previous_image, rtol=1e-5, atol=1e-8)

        previous_flux, previous_image = flux, image

def test_tempgrad_asymmetry_and_rotation():
    """Test skew and rotation with oimTempGrad in L and N bands."""

    base = dict(
        rin=0.15,
        rout=1.5,
        r0=1.0,
        T0=600.0,
        q=-0.5,
        p=-0.5,
        log_sigma0=-3.0,
        compute_sigma0=False,
        kappa_abs=300.0,
        dist=154.1,
        dim=256,
        elong=1.25,
        pa=108.0,
    )

    rotating = oim.oimTempGrad(
        **base,
        asymmetric=True,
        modulation=1,
        skw1=0.35,
        skwPa1=30.0,
        rotPeriod=20.0,
        rotT0=60000.0,
    )

    # Five days gives a 90 degree turn
    static = oim.oimTempGrad(
        **base,
        asymmetric=True,
        modulation=1,
        skw1=0.35,
        skwPa1=120.0,
    )

    symmetric = oim.oimTempGrad(
        **base,
        asymmetric=False,
    )

    u = np.array([0.0, 1e7, 2e7, -3e7])
    v = np.array([0.0, 2e7, -1e7, 2e7])

    epoch = 60005.0

    # L and N bands
    for wavelength in [3.5e-6, 10.0e-6]:

        wl = np.full(u.size, wavelength)
        t = np.full(u.size, epoch)

        f_rot = rotating.getComplexCoherentFlux(u, v, wl, t)
        f_static = static.getComplexCoherentFlux(u, v, wl, t)
        f_sym = symmetric.getComplexCoherentFlux(u, v, wl, t)

        assert np.all(np.isfinite(f_rot))
        assert np.abs(f_rot[0]) > 0

        # Same phase, same result
        np.testing.assert_allclose(
            f_rot,
            f_static,
            rtol=1e-11,
            atol=1e-10,
        )

        # Same total flux
        np.testing.assert_allclose(
            f_rot[0],
            f_sym[0],
            rtol=1e-12,
            atol=0,
        )

        # Skew changes the visibility
        v_rot = f_rot / f_rot[0]
        v_sym = f_sym / f_sym[0]

        assert np.max(np.abs(v_rot[1:] - v_sym[1:])) > 1e-4

        # Compare images
        def image(component):
            return component.getImage(
                dim=161,
                pixSize=0.16,
                wl=np.array([wavelength]),
                t=np.array([epoch]),
            )[0, 0]

        im_rot = image(rotating)
        im_static = image(static)
        im_sym = image(symmetric)

        assert np.all(np.isfinite(im_rot))
        assert np.all(np.isfinite(im_static))
        assert np.all(np.isfinite(im_sym))

        assert np.all(np.isfinite(f_rot))
        assert np.all(np.isfinite(f_static))
        assert np.all(np.isfinite(f_sym))

        np.testing.assert_allclose(
            im_rot,
            im_static,
            rtol=1e-12,
            atol=1e-12,
        )

        assert np.max(np.abs(im_rot - im_sym)) > 1e-8
