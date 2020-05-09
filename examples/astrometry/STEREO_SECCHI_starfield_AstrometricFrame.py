"""
==================================================================================
Identifying stars in a STEREO/SECCHI COR2 coronagraph image using AstrometricFrame
==================================================================================

Since the field of view in SECCHI COR2 images can span 2 to 15 solar radii,
we often observe stars in these image data. In this example, we will use the
`Astroquery package <https://astroquery.readthedocs.io/en/latest/>`__ to query the
`VizieR star catalog <http://vizier.u-strasbg.fr/viz-bin/VizieR>`__ for stars observed
by the `Gaia satellite <https://sci.esa.int/web/gaia/>`__ within the SECCHI COR2 field of view.
Then we will use the coordinates framework in SunPy and AstroPy to transform the coordinates
returned by VizieR into SECCHI COR2 image coordinates. As a bonus, we'll also identify Mars."

This requires the installation of the `astroquery <https://astroquery.readthedocs.io/en/latest/>`__
package, which can be installed on top of the existing sunpy conda
environment: ``conda install -c astropy astroquery`` and an active internet connection.
"""
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier

import astropy.units as u
from astropy.coordinates import Distance, SkyCoord

import sunpy.map
from sunpy.coordinates import AstrometricFrame, get_body_heliographic_stonyhurst
from sunpy.net import helioviewer

###############################################################################
# Let's download a STEREO-A SECCHI COR2 image from Helioviewer which provide
# pre-processed images and load it into a Map.

hv = helioviewer.HelioviewerClient()
f = hv.download_jp2('2014/05/15 07:54', observatory='STEREO_A',
                    instrument='SECCHI', detector='COR2')
cor2 = sunpy.map.Map(f)

###############################################################################
# To efficiently search the star field, we need to know what stars are near the
# Sun as observed by STEREO. We need the vector that points from STEREO to the Sun.
# The location of STEREO in HCRS provides the Sun-to-STEREO vector.

sun_to_stereo = cor2.observer_coordinate.transform_to('hcrs')

###############################################################################
# We next reflect the vector to get our search vector which points from STEREO
# to the Sun.

stereo_to_sun = SkyCoord(-sun_to_stereo.spherical, obstime=sun_to_stereo.obstime, frame='hcrs')

###############################################################################
# Let's look up bright stars using the Vizier search capability provided by
# astroquery.
# We will search the GAIA2 star catalog for stars with magnitude
# brighter than 7.

vv = Vizier(columns=['**'], row_limit=-1, column_filters={'Gmag': '<7'}, timeout=1200)
vv.ROW_LIMIT = -1
result = vv.query_region(stereo_to_sun, radius=4 * u.deg, catalog='I/345/gaia2')

###############################################################################
# Let's see how many stars we've found.

print(len(result[0]))

###############################################################################
# Now we load all stars into an array coordinate and transform it into the COR2
# image coordinates.  The reference epoch for the star positions is J2015.5,
# and we update these positions to the date of the COR2 observation.

tbl_crds = SkyCoord(ra=result[0]['RA_ICRS'],
                    dec=result[0]['DE_ICRS'],
                    distance=Distance(parallax=u.Quantity(result[0]['Plx'])),
                    pm_ra_cosdec=result[0]['pmRA'],
                    pm_dec=result[0]['pmDE'],
                    radial_velocity=result[0]['RV'],
                    frame='icrs', obstime='J2015.5')
cor2_ast_frame = AstrometricFrame(base=cor2.coordinate_frame,
                                  observer=cor2.observer_coordinate,
                                  ref_epoch=tbl_crds.obstime)
hpc_coords = tbl_crds.transform_to(cor2_ast_frame)

###############################################################################
# One of the bright features is actually Mars, so let's also get that coordinate.

mars = get_body_heliographic_stonyhurst('mars', cor2.date)
mars_hpc = mars.transform_to(cor2_ast_frame)

###############################################################################
# Let's plot the results.

ax = plt.subplot(projection=cor2)

# Let's tweak the axis to show in degrees instead of arcsec
lon, lat = ax.coords
lon.set_major_formatter('d.dd')
lat.set_major_formatter('d.dd')

cor2.plot(axes=ax, vmin=0, vmax=600)
cor2.draw_limb()

# Plot the position of Mars
ax.plot_coord(mars_hpc.as_base(), 's', color='white', fillstyle='none', markersize=12, label='Mars')
# Plot all of the stars
ax.plot_coord(hpc_coords.as_base(), 'o', color='white', fillstyle='none')
plt.legend()
plt.show()