=====
Usage
=====

To use Xspharm in a project::

    import xspharm
    xsp = xspharm.xspharm(grid_ds, gridtype='regular')
    
    # tri_truncate field using total wavenumber ntrunc
    field_t = xsp.truncate(field_ds, ntrunc=24)

    # smoothing (tapering) the fields
    field_filter = xsp.exp_taper(field_ds, ntrunc=24, r=2)

    

