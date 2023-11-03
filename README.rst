=================
 Xspharm Package
=================

.. image:: https://img.shields.io/pypi/v/xspharm.svg
    :target: https://pypi.python.org/pypi/xspharm
    :alt: PyPI Version

.. image:: https://readthedocs.org/projects/xspharm/badge/?version=latest
    :target: https://xspharm.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Overview
--------

Xspharm is an xarray-compatible library that facilitates spherical harmonic transforms. It leverages the computational efficiency of `pyspharm` and the convenience of `xarray` data structures to provide an intuitive interface for processing geospatial data on a spherical domain.

* **License**: Distributed under the BSD License.
* **Documentation**: Comprehensive documentation available at https://xspharm.readthedocs.io.

Features
--------

Xspharm provides a suite of methods to manipulate and transform geospatial datasets:

- `truncate`: Reduces the resolution of a data variable or an entire dataset to a specified spherical harmonic wavenumber.
- `exp_taper`: Applies tapering to mitigate the Gibbs phenomenon in spherical harmonic coefficients.
- `uv2sfvp`: Transforms zonal (`u`) and meridional (`v`) wind components into streamfunction (`sf`) and velocity potential (`vp`).
- `uv2vordiv`: Converts zonal (`u`) and meridional (`v`) wind components to vorticity and divergence fields.
- `uv2absvor`: Changes zonal (`u`) and meridional (`v`) wind components to absolute vorticity.
- `sf2uv`: Derives rotational wind components from a given streamfunction.
- `vp2uv`: Obtains divergent wind components from velocity potential.
- `sfvp2uv`: Integrates streamfunction and velocity potential to produce zonal (`u`) and meridional (`v`) wind components.

Credits
-------

Special thanks to `Jeff Whitaker`_ for sharing the `pyspharm`_ library, which forms the foundation of this package's capabilities; and to `Andrew Dawson`_ for inspiring with the `windspharm`_ library.

.. _Jeff Whitaker: https://github.com/jswhit
.. _Andrew Dawson: https://github.com/ajdawson
.. _pyspharm: https://github.com/jswhit/pyspharm
.. _windspharm: https://github.com/ajdawson/windspharm

