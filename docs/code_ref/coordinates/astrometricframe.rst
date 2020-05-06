.. _sunpy-coordinates-astrometricframe:

Astrometric frames
==================

Examples for `~sunpy.coordinates.metaframes.AstrometricFrame`

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
  >>> venus.transform_to(hgs_ast_frame)
  <AstrometricHeliographicStonyhurst Coordinate (base=<HeliographicStonyhurst Frame (obstime=2012-06-06T04:07:29.000)>, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>): (lon, lat, radius) in (deg, deg, AU)
      (0.07084923, 0.0520573, 0.72605477)
   (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
      (0.01084266, -0.0025724, -0.07821503)>

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
      (7.835757e-15, -0.00766698, 1.01475668)>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (641.43738743, 540.71510388, 0.28870525)
   (d_Tx, d_Ty, d_distance) in (arcsec / s, arcsec / s, km / s)
      (0.12666774, -0.00647582, 0.1918771)>
  >>> venus_hpc_ast.Tx, venus_hpc_ast.Ty
  (<Longitude 641.43738743 arcsec>, <Latitude 540.71510388 arcsec>)

Compare against the light-travel-time functionality of `~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`::

  >>> get_body_heliographic_stonyhurst('venus', obstime, observer=earth).transform_to(hpc_frame)
  INFO: Apparent body location accounts for 144.07 seconds of light travel time [sunpy.coordinates.ephemeris]
  <Helioprojective Coordinate (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (641.43758224, 540.71512132, 0.28870525)>

Transforming back to the base frame will just undo the astrometric adjustment::

  >>> venus_hpc_ast.transform_to(hpc_frame)
  <Helioprojective Coordinate (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (665.39449753, 542.33097368, 0.28870521)
   (d_Tx, d_Ty, d_distance) in (arcsec / s, arcsec / s, km / s)
      (0.12666765, -0.00647588, 0.19590831)>

Instead, casting to a real (non-astrometric) coordinate frame::

  >>> venus_hpc_ast.as_base()
  <Helioprojective Coordinate (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
      (7.835757e-15, -0.00766698, 1.01475668)>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (641.43738743, 540.71510388, 0.28870525)
   (d_Tx, d_Ty, d_distance) in (arcsec / s, arcsec / s, km / s)
      (0.12666774, -0.00647582, 0.1918771)>

