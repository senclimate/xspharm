import numpy as np
import xarray as xr
from spharm import Spharmt

class xspharm:
    """
    Xarray interface to spherical harmonic transforms using pyspharm.

    Attributes:
        nlat (int): Number of latitude points.
        nlon (int): Number of longitude points.
        fvor (xarray DataArray) : planetary vorticity 
        s (Spharmt): Spharmt object for spherical harmonic operations.

    Methods:
        truncate: Truncate a data variable or entire dataset to a specific wavenumber.
        exp_taper: Apply tapering to spherical harmonic coefficients.
        uv2sfvp: Convert zonal and meridional wind components to streamfunction and velocity potential.
        uv2vordiv: Convert zonal and meridional wind components to vorticity and divergence.
        uv2absvor: Convert zonal and meridional wind components to absolute vorticity.
        sf2uv: Convert streamfunction to rotational wind components.
        vp2uv: Convert velocity potential to divergent wind components.
        sfvp2uv: Convert streamfunction and velocity potential to zonal and meridional wind components.
    """

    def __init__(self, grid_ds,
                 gridtype='regular',
                 rsphere=6.3712e6, 
                 omega=7.292e-5, 
                 legfunc='stored'):
        """
        Initialize the XSpharm class with the given dataset and spherical harmonic transform parameters.

        Args:
            grid_ds (xarray.Dataset): Dataset containing the grid information.
            gridtype (str): Type of grid, either 'regular' or 'gaussian'.
            rsphere (float): Radius of the sphere, defaults to Earth's radius in meters.
            omega (float): Rotation rate of the planet, defaults to Earth's rotation rate.
            legfunc (str): Legendre function computation method, either 'stored' or 'computed'.
        """
        self.nlat = len(grid_ds['lat'])
        self.nlon = len(grid_ds['lon'])
        self.fvor = 2. * omega * np.sin(np.deg2rad(grid_ds['lat'])) 
        self.s = Spharmt(self.nlon, self.nlat, gridtype=gridtype, rsphere=rsphere, legfunc=legfunc)

    def truncate(self, ds_or_var, ntrunc=None):
        """
        Truncate a data variable or entire dataset.

        Args:
        - ds_or_var (xr.DataArray or xr.Dataset): Input data.
        - ntrunc (int, optional): Truncation wavenumber. Default is None.

        Returns:
        xr.DataArray or xr.Dataset: Truncated data.
        """
        if isinstance(ds_or_var, xr.Dataset):
            # Create a new dataset by applying truncation to each data variable
            truncated_data = {}
            for var_name, var_data in ds_or_var.data_vars.items():
                truncated_data[var_name] = self._truncate_var(var_data, ntrunc)
            return xr.Dataset(truncated_data)
        elif isinstance(ds_or_var, xr.DataArray):
            return self._truncate_var(ds_or_var, ntrunc)
        else:
            raise ValueError("Input should be an xarray Dataset or DataArray.")

    def _truncate_var(self, var_ds, ntrunc=None):
        other_dims = _get_other_dims(var_ds)
        var_data = _transpose_to_spharm(var_ds)
        spec_data = self.s.grdtospec(var_data.values, ntrunc=ntrunc)
        grid_data = var_data.copy(data=self.s.spectogrd(spec_data))
        return _transpose_from_spharm(grid_data, other_dims)

    def exp_taper(self, ds_or_var, ntrunc=None, r=2):
        """
        Taper (filter) the spherical harmonic coefficients using the equation provided.

        Args:
        - ds_or_var (xr.DataArray or xr.Dataset): Input data.
        - ntrunc (int): Truncation wavenumber N.
        - r (float): User-defined parameter.

        Returns:
        xr.DataArray or xr.Dataset: Tapered data.
        """
        if isinstance(ds_or_var, xr.Dataset):
            tapered_data = {}
            for var_name, var_data in ds_or_var.data_vars.items():
                tapered_data[var_name] = self._exp_taper_var(var_data, ntrunc, r)
            return xr.Dataset(tapered_data)
        elif isinstance(ds_or_var, xr.DataArray):
            return self._exp_taper_var(ds_or_var, ntrunc, r)
        else:
            raise ValueError("Input should be an xarray Dataset or DataArray.")

    def _exp_taper_var(self, var_ds, ntrunc, r):
        other_dims = _get_other_dims(var_ds)
        var_data = _transpose_to_spharm(var_ds)

        # Compute spherical harmonic coefficients
        spec_data = self.s.grdtospec(var_data.values)

        # Compute l and m values based on triangular layout
        l_values, m_values = [], []
        for m in range(self.nlat-1 + 1):
            for l in range(m, self.nlat-1 + 1):
                l_values.append(l)
                m_values.append(m)

        l_values = np.array(l_values, dtype=float)
        m_values = np.array(m_values, dtype=float)
        total_wavenumber = np.sqrt(l_values * (l_values + 1) + m_values**2)

        # Apply tapering based on the given equation
        taper_filter = np.exp(-((total_wavenumber*(total_wavenumber+1))/(ntrunc*(ntrunc+1)))**r)
        spec_data *= taper_filter[:, np.newaxis] if spec_data.ndim == 2 else taper_filter

        # Convert back to grid data
        grid_data = var_data.copy(data=self.s.spectogrd(spec_data))
        return _transpose_from_spharm(grid_data, other_dims)

    def uv2sfvp(self, u_ds, v_ds, ntrunc=None):
        """Streamfunction and velocity potential.
        
        Inputs:
        u, v: zonal and meridional winds
        ntrunc: Truncation limit (triangular truncation) for the spherical harmonic computation.

        Returns:
        sf, vp
            The streamfunction and velocity potential respectively.
        """
        other_dims = _get_other_dims(u_ds)
        u_data = _transpose_to_spharm(u_ds)
        v_data = _transpose_to_spharm(v_ds)
        psi_data, chi_data = self.s.getpsichi(u_data.values, v_data.values, ntrunc=ntrunc)

        psi_ds = _transpose_from_spharm(u_data.copy(data=psi_data), other_dims)
        chi_ds = _transpose_from_spharm(u_data.copy(data=chi_data), other_dims)
        psi_ds.attrs['long_name'] = 'streamfunction'
        psi_ds.attrs['units'] = 'm**2/s'
        chi_ds.attrs['long_name'] = 'velocity potential'
        chi_ds.attrs['units'] = 'm**2/s'
        return xr.Dataset({'sf': psi_ds, 'vp': chi_ds})

    def uv2vordiv(self, u_ds, v_ds, ntrunc=None):
        """Computes the vorticity and divergence via spherical harmonics, given the u and v wind components
        
        Inputs:
        u, v: zonal and meridional winds
        ntrunc: Truncation limit (triangular truncation) for the spherical harmonic computation.

        Returns:
        vor, div
            The vorticity and divergence respectively.
        """
        other_dims = _get_other_dims(u_ds)
        u_data = _transpose_to_spharm(u_ds)
        v_data = _transpose_to_spharm(v_ds)

        vor_spec, div_spec = self.s.getvrtdivspec(u_data.values, v_data.values, ntrunc=ntrunc)
        vor_ds = _transpose_from_spharm( u_data.copy(data=self.s.spectogrd(vor_spec)), other_dims)
        div_ds = _transpose_from_spharm( u_data.copy(data=self.s.spectogrd(div_spec)), other_dims)
        vor_ds.attrs['long_name'] = 'Vorticity'
        vor_ds.attrs['units'] = '1/s'
        div_ds.attrs['long_name'] = 'Divergence'
        div_ds.attrs['units'] = '1/s'
        return xr.Dataset({'vor': vor_ds, 'div': div_ds})

    def uv2absvor(self, u_ds, v_ds, ntrunc=None):
        """
        Computes the absolute vorticity via spherical harmonics, given the u and v wind components
        
        Inputs:
        u, v: zonal and meridional winds
        ntrunc: Truncation limit (triangular truncation) for the spherical harmonic computation.

        Returns:
        vor, div
            The vorticity and divergence respectively.
        """
        other_dims = _get_other_dims(u_ds)
        u_data = _transpose_to_spharm(u_ds)
        v_data = _transpose_to_spharm(v_ds)

        vor_spec, _ = self.s.getvrtdivspec(u_data.values, v_data.values, ntrunc=ntrunc)
        vor_ds = _transpose_from_spharm( u_data.copy(data=self.s.spectogrd(vor_spec)), other_dims)
        
        # add planetary vorticity (Coriolis effect) 
        vor_ds = vor_ds + self.fvor
        vor_ds.attrs['long_name'] = 'absolute vorticity'
        vor_ds.attrs['units'] = '1/s'
        return vor_ds

    def sf2uv(self, sf_ds, ntrunc=None):
        """Computes the rotational wind components via spherical harmonics, given an array containing streamfunction
        
        Inputs:
        sf_ds streamfunction
        ntrunc: Truncation limit (triangular truncation) for the spherical harmonic computation.

        Returns:
            u, v: zonal and meridional winds
        """
        other_dims = _get_other_dims(sf_ds)
        sf_data = _transpose_to_spharm(sf_ds)
        psi_spec = self.s.grdtospec(sf_data, ntrunc=ntrunc)
        vpsi_data, upsi_data = self.s.getgrad(psi_spec)
        
        u_ds = _transpose_from_spharm( sf_data.copy(data=-upsi_data), other_dims)
        v_ds = _transpose_from_spharm( sf_data.copy(data=vpsi_data), other_dims)
        u_ds.attrs['long_name'] = 'rotational component of U wind'
        u_ds.attrs['units'] = 'm/s'
        v_ds.attrs['long_name'] = 'rotational component of V wind'
        v_ds.attrs['units'] = 'm/s'
        return xr.Dataset({'u_rot': u_ds, 'v_rot': v_ds})

    def vp2uv(self, vp_ds, ntrunc=None):
        """Computes the divergent wind components via spherical harmonics, given an array containing velocity potential 
        
        Inputs:
        vp_ds velocity potential.
        ntrunc: Truncation limit (triangular truncation) for the spherical harmonic computation.

        Returns:
            u, v: zonal and meridional winds
        """
        other_dims = _get_other_dims(vp_ds)
        vp_data = _transpose_to_spharm(vp_ds)
        chi_spec = self.s.grdtospec(vp_data, ntrunc=ntrunc)
        udiv_data, vdiv_data = self.s.getgrad(chi_spec)
        
        u_ds = _transpose_from_spharm( vp_data.copy(data=udiv_data), other_dims)
        v_ds = _transpose_from_spharm( vp_data.copy(data=vdiv_data), other_dims)
        u_ds.attrs['long_name'] = 'divergent component of U wind'
        u_ds.attrs['units'] = 'm/s'
        v_ds.attrs['long_name'] = 'divergent component of V wind'
        v_ds.attrs['units'] = 'm/s'
        return xr.Dataset({'u_div': u_ds, 'v_div': v_ds})

    def sfvp2uv(self, sf_ds, vp_ds, ntrunc=None):
        """
        Computes the wind components via spherical harmonics, given streamfunction and velocity potential.
        
        Inputs:
        sf_ds: streamfunction
        vp_ds: velocity potential
        ntrunc: Truncation limit (triangular truncation) for the spherical harmonic computation.
        
        Returns:
            u, v: zonal and meridional winds
        """
        other_dims = _get_other_dims(sf_ds)
        sf_data = _transpose_to_spharm(sf_ds)
        vp_data = _transpose_to_spharm(vp_ds)

        psi_spec = self.s.grdtospec(sf_data, ntrunc=ntrunc)
        vpsi_data, upsi_data = self.s.getgrad(psi_spec)
        chi_spec = self.s.grdtospec(vp_data, ntrunc=ntrunc)
        udiv_data, vdiv_data = self.s.getgrad(chi_spec)

        u_ds = _transpose_from_spharm( sf_data.copy(data=udiv_data-upsi_data), other_dims)
        v_ds = _transpose_from_spharm( sf_data.copy(data=vdiv_data+vpsi_data), other_dims)
        u_ds.attrs['long_name'] = 'U wind'
        u_ds.attrs['units'] = 'm/s'
        v_ds.attrs['long_name'] = 'V wind'
        v_ds.attrs['units'] = 'm/s'
        return xr.Dataset({'u': u_ds, 'v': v_ds})


def _get_other_dims(input_data):
    '''
    Get other dimensions besides 'lat' and 'lon'.
    '''
    spatial_dims = ['lat', 'lon']
    other_dims = [dim for dim in input_data.dims if dim not in spatial_dims]
    return other_dims


def _transpose_to_spharm(input_data):
    """
    Transpose data to fit spharm's expected layout: [nlat, nlon, nt].
    """
    other_dims = _get_other_dims(input_data)

    if len(other_dims) == 0:
        return input_data
    elif len(other_dims) == 1:
        return input_data.transpose('lat', 'lon', other_dims[0])
    elif len(other_dims) >= 2:
        return input_data.stack(nt=other_dims).transpose('lat', 'lon', 'nt')
    else:
        raise ValueError("Input data dimensions not supported.")


def _transpose_from_spharm(sp_data, other_dims):
    """
    Transpose data back from spharm's layout.
    """
    if len(other_dims) == 0:
        return sp_data
    elif len(other_dims) == 1:
        return sp_data.transpose(other_dims[0], 'lat', 'lon')
    elif len(other_dims) >= 2:
        return sp_data.unstack('nt').transpose(*other_dims, 'lat', 'lon')
    else:
        raise ValueError("Output data dimensions not supported.")