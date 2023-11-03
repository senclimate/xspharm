=======
Xspharm
=======


.. image:: https://img.shields.io/pypi/v/xspharm.svg
        :target: https://pypi.python.org/pypi/xspharm

.. image:: https://readthedocs.org/projects/xspharm/badge/?version=latest
        :target: https://xspharm.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




xarray interface to spherical harmonic transform


* Free software: BSD license
* Documentation: https://xspharm.readthedocs.io.


Features
--------

* Methods
truncate: Truncate a data variable or entire dataset to a specific wavenumber.
exp_taper: Apply tapering to spherical harmonic coefficients.
uv2sfvp: Convert zonal and meridional wind components to streamfunction and velocity potential.
uv2vordiv: Convert zonal and meridional wind components to vorticity and divergence.
uv2absvor: Convert zonal and meridional wind components to absolute vorticity.
sf2uv: Convert streamfunction to rotational wind components.
vp2uv: Convert velocity potential to divergent wind components.
sfvp2uv: Convert streamfunction and velocity potential to zonal and meridional wind components.


Credits
-------

Thanks go to _Jeff Whitaker and _Andrew Dawson for sharing _pyspharm and _windspharm libraries, respectively.

.. _Jeff Whitaker: https://github.com/jswhit
.. _pyspharm: https://github.com/jswhit/pyspharm
.. _Andrew Dawson: https://github.com/ajdawson
.. _windspharm: https://github.com/ajdawson/windspharm
