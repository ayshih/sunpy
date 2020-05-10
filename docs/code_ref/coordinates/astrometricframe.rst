.. _sunpy-coordinates-astrometricframe:

Astrometric frames
==================

Tutorial for `~sunpy.coordinates.metaframes.AstrometricFrame`

Astrometric location of a solar-system body
-------------------------------------------

`~sunpy.coordinates.metaframes.AstrometricFrame` computes the astrometric location of a solar-system body by assuming that it is moving in a straight line with its instantaneous velocity.

Instantaneous location of Venus::

  >>> from sunpy.coordinates import (AstrometricFrame, HeliographicStonyhurst,
  ...                                Helioprojective, get_body_heliographic_stonyhurst)
  >>> obstime = '2012-06-06 04:07:29'
  >>> earth = get_body_heliographic_stonyhurst('earth', obstime)
  >>> venus = get_body_heliographic_stonyhurst('venus', obstime, include_velocity=True)
  >>> venus
  <HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (0.07349535, 0.05223575, 0.72605496)
   (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
      (0.01084266, -0.0025724, -0.07795555)>

Astrometric location of Venus::

  >>> hgs_ast_frame = AstrometricFrame(base=HeliographicStonyhurst(obstime=obstime), observer=earth)
  >>> venus_ast = venus.transform_to(hgs_ast_frame)
  >>> venus_ast
  <AstrometricHeliographicStonyhurst Coordinate (base=<HeliographicStonyhurst Frame (obstime=2012-06-06T04:07:29.000)>, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>, ref_epoch=J2000.000): (lon, lat, radius) in (deg, deg, AU)
      (0.07084923, 0.0520573, 0.72605477)
   (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
      (0.01084266, -0.0025724, -0.07821503)>

Astrometric location of a cosmic object
---------------------------------------

`~sunpy.coordinates.metaframes.AstrometricFrame` internally calls :meth:`astropy.coordinates.SkyCoord.apply_space_motion` to propagate the position of cosmic objects.
Example::

  >>> import astropy.units as u
  >>> from astropy.time import Time
  >>> from astropy.coordinates import SkyCoord
  >>> c = SkyCoord(ra=10*u.degree, dec=45*u.degree, distance=100*u.pc,
  ...              pm_ra_cosdec=200*u.mas/u.yr, pm_dec=-300*u.mas/u.yr,
  ...              radial_velocity=100*u.km/u.s,
  ...              frame='icrs', obstime=Time('J2010'))

  >>> new_obstime = Time('2011-12-13 14:15:16')
  >>> hgs_frame = HeliographicStonyhurst(obstime=new_obstime)
  >>> c.transform_to(hgs_frame)
  <SkyCoord (HeliographicStonyhurst: obstime=2011-12-13 14:15:16.000): (lon, lat, radius) in (deg, deg, pc)
      (-47.22670834, 41.92004583, 100.00000001)
   (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
      (-2.44709473e-09, -1.2216929e-08, 119.03469026)>

  >>> cc = c.apply_space_motion(new_obstime=new_obstime)
  >>> cc.transform_to(hgs_frame)
  <SkyCoord (HeliographicStonyhurst: obstime=2011-12-13 14:15:16.000): (lon, lat, radius) in (deg, deg, pc)
      (-47.22671914, 41.91985087, 100.00019913)
   (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
      (-2.44707881e-09, -1.22168778e-08, 119.03531576)>

Using `~sunpy.coordinates.metaframes.AstrometricFrame` instead, and the observer does not need to be specified::

  >>> ccc = c.transform_to(AstrometricFrame(base=hgs_frame, ref_epoch='J2010'))
  >>> ccc
  <SkyCoord (AstrometricHeliographicStonyhurst: base=<HeliographicStonyhurst Frame (obstime=2011-12-13 14:15:16.000)>, observer=None, ref_epoch=J2010.000): (lon, lat, radius) in (deg, deg, pc)
      (-47.22671904, 41.91985104, 100.0001998)
   (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
      (-2.44149313e-09, -1.22067667e-08, 119.36867142)>

The discrepancies compared to the earlier approach are due to the size of the finite-difference time step used to calculate the velocity (see also :ref:`astropy-coordinate-finite-difference-velocities`).
The default time step is 1 second, but the motion of a star is very small compared to its distance on this short of a time scale.
If the time step is changed to 1 year, the discrepancy becomes negligible.
We use the context manager :func:`~sunpy.coordinates.transformations.impose_finite_difference_dt` to change the finite-difference time step (``finite_difference_dt``) for all transformations::

  >>> from sunpy.coordinates import impose_finite_difference_dt
  >>> with impose_finite_difference_dt(1*u.yr):
  ...     cccc = c.transform_to(AstrometricFrame(base=hgs_frame, ref_epoch='J2010'))
  >>> cccc
  <SkyCoord (AstrometricHeliographicStonyhurst: base=<HeliographicStonyhurst Frame (obstime=2011-12-13 14:15:16.000)>, observer=None, ref_epoch=J2010.000): (lon, lat, radius) in (deg, deg, pc)
      (-47.22671914, 41.91985087, 100.00019913)
   (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
      (-2.44707881e-09, -1.22168778e-08, 119.03531577)>

Transforming back reverses the propagation (again, using the context manager)::

  >>> with impose_finite_difference_dt(1*u.yr):
  ...    print(cccc.transform_to(hgs_frame))
  <SkyCoord (HeliographicStonyhurst: obstime=2011-12-13 14:15:16.000): (lon, lat, radius) in (deg, deg, pc)
      (-47.22670834, 41.92004583, 100.00000001)
   (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
      (-2.44709473e-09, -1.2216929e-08, 119.03469026)>

Astrometric location in an observer-based frame
-----------------------------------------------

Instantaneous location of Venus in `~sunpy.coordinates.frames.Helioprojective`::

  >>> hpc_frame = Helioprojective(observer=earth, obstime=obstime)
  >>> venus_hpc = venus.transform_to(hpc_frame)
  >>> venus_hpc
  <Helioprojective Coordinate (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (665.39449468, 542.33097349, 0.28870521)
   (d_Tx, d_Ty, d_distance) in (arcsec / s, arcsec / s, km / s)
      (0.12666765, -0.00647588, 0.19590823)>

  >>> venus_hpc.Tx, venus_hpc.Ty
  (<Longitude 665.39449468 arcsec>, <Latitude 542.33097349 arcsec>)

Astrometric location of Venus in `~sunpy.coordinates.frames.Helioprojective`::

  >>> hpc_ast_frame = AstrometricFrame(base=hpc_frame)  # The observer is pulled from `hpc_frame`
  >>> venus_hpc_ast = venus_hpc.transform_to(hpc_ast_frame)
  >>> venus_hpc_ast
  <AstrometricHelioprojective Coordinate (base=<Helioprojective Frame (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>)>, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>, ref_epoch=J2000.000): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (641.43738458, 540.71510369, 0.28870525)
   (d_Tx, d_Ty, d_distance) in (arcsec / s, arcsec / s, km / s)
      (0.12666774, -0.00647582, 0.1918771)>

  >>> venus_hpc_ast.Tx, venus_hpc_ast.Ty
  (<Longitude 641.43738458 arcsec>, <Latitude 540.71510369 arcsec>)

Obtaining real coordinates
--------------------------

Transforming back to the base frame will just undo the astrometric adjustment::

  >>> venus_hpc_ast.transform_to(hpc_frame)
  <Helioprojective Coordinate (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (665.39449468, 542.33097349, 0.28870521)
   (d_Tx, d_Ty, d_distance) in (arcsec / s, arcsec / s, km / s)
      (0.12666765, -0.00647588, 0.19590826)>

Instead, casting to a real (non-astrometric) coordinate frame::

  >>> venus_hpc_ast.as_base()
  <Helioprojective Coordinate (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (641.43738458, 540.71510369, 0.28870525)
   (d_Tx, d_Ty, d_distance) in (arcsec / s, arcsec / s, km / s)
      (0.12666774, -0.00647582, 0.1918771)>

Comparison of accuracy
----------------------

Compare against the light-travel-time functionality of `~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`::

  >>> venus_accurate = SkyCoord(get_body_heliographic_stonyhurst('venus', obstime, observer=earth))
  INFO: Apparent body location accounts for 144.07 seconds of light travel time [sunpy.coordinates.ephemeris]
  >>> SkyCoord(venus_ast).separation_3d(venus_accurate).to('km')
  <Distance 0.16215772 km>

For this coordinate, which is during a transit of Venus across the Sun, this distance discrepancy is primarily along the line of sight, so the angular discrepancy in the sky is very small::

  >>> venus_hpc_accurate = venus_accurate.transform_to(hpc_frame)
  >>> SkyCoord(venus_hpc_ast.as_base()).separation(venus_hpc_accurate).to('arcsec')
  <Angle 0.00019845 arcsec>

Example uses of AstrometricFrame
--------------------------------

Here are the examples in our gallery that use `~sunpy.coordinates.metaframes.AstrometricFrame`:

.. include:: ../../generated/modules/sunpy.coordinates.AstrometricFrame.examples
    :start-line: 5

.. raw:: html

    <div class="sphx-glr-clear"></div>
