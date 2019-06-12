# -*- coding: utf-8 -*-
"""
Coordinate Transformation Functions

This module contains the functions for converting one
`sunpy.coordinates.frames` object to another.

.. warning::

  The functions in this submodule should never be called directly, transforming
  between coordinate frames should be done using the ``.transform_to`` methods
  on `~astropy.coordinates.BaseCoordinateFrame` or
  `~astropy.coordinates.SkyCoord` instances.

"""
from copy import deepcopy

import numpy as np

import astropy.units as u
from astropy.coordinates import ICRS, HCRS, ConvertError, BaseCoordinateFrame, get_body_barycentric
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.representation import (CartesianRepresentation, SphericalRepresentation,
                                                UnitSphericalRepresentation)
from astropy.coordinates.transformations import (FunctionTransform, DynamicMatrixTransform,
                                                 FunctionTransformWithFiniteDifference)
from astropy.coordinates.matrix_utilities import matrix_product, rotation_matrix, matrix_transpose
# Versions of Astropy that do not have HeliocentricMeanEcliptic have the same frame
# with the incorrect name HeliocentricTrueEcliptic
try:
    from astropy.coordinates import HeliocentricMeanEcliptic
except ImportError:
    from astropy.coordinates import HeliocentricTrueEcliptic as HeliocentricMeanEcliptic

from sunpy.sun import constants

from .frames import (Heliocentric, Helioprojective, HeliographicCarrington, HeliographicStonyhurst,
                     HeliocentricEarthEcliptic, GeocentricSolarEcliptic, HeliocentricInertial)

try:
    from astropy.coordinates.builtin_frames import _make_transform_graph_docs as make_transform_graph_docs
except ImportError:
    from astropy.coordinates import make_transform_graph_docs as _make_transform_graph_docs
    make_transform_graph_docs = lambda: _make_transform_graph_docs(frame_transform_graph)


RSUN_METERS = constants.get('radius').si.to(u.m)

__all__ = ['hgs_to_hgc', 'hgc_to_hgs', 'hcc_to_hpc',
           'hpc_to_hcc', 'hcc_to_hgs', 'hgs_to_hcc',
           'hpc_to_hpc',
           'hcrs_to_hgs', 'hgs_to_hcrs',
           'hgs_to_hgs', 'hgc_to_hgc', 'hcc_to_hcc',
           'hme_to_hee', 'hee_to_hme', 'hee_to_hee',
           'hee_to_gse', 'gse_to_hee', 'gse_to_gse',
           'hme_to_hci', 'hci_to_hme', 'hci_to_hci']


def _carrington_offset(obstime):
    """
    Calculate the HG Longitude offest based on a time
    """
    if obstime is None:
        raise ValueError("To perform this transformation the coordinate"
                         " Frame needs a obstime Attribute")

    # Import here to avoid a circular import
    from .sun import L0
    return L0(obstime)


def _observers_are_equal(obs_1, obs_2, string_ok=False):
    if string_ok:
        if obs_1 == obs_2:
            return True
    if not (isinstance(obs_1, BaseCoordinateFrame) and isinstance(obs_2, BaseCoordinateFrame)):
        raise ValueError("To compare two observers, both must be instances of BaseCoordinateFrame. "
                         "Cannot compare two observers {} and {}.".format(obs_1, obs_2))
    return (u.allclose(obs_1.lat, obs_2.lat) and
            u.allclose(obs_1.lon, obs_2.lon) and
            u.allclose(obs_1.radius, obs_2.radius))


# =============================================================================
# ------------------------- Transformation Framework --------------------------
# =============================================================================


@frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst,
                                 HeliographicCarrington)
def hgs_to_hgc(hgscoord, hgcframe):
    """
    Transform from Heliographic Stonyhurst to Heliograpic Carrington.
    """
    if hgcframe.obstime is None or np.any(hgcframe.obstime != hgscoord.obstime):
        raise ValueError("Can not transform from Heliographic Stonyhurst to "
                         "Heliographic Carrington, unless both frames have matching obstime.")

    c_lon = hgscoord.spherical.lon + _carrington_offset(hgscoord.obstime).to(u.deg)
    representation = SphericalRepresentation(c_lon, hgscoord.spherical.lat,
                                             hgscoord.spherical.distance)
    hgcframe = hgcframe.__class__(obstime=hgscoord.obstime)

    return hgcframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HeliographicCarrington,
                                 HeliographicStonyhurst)
def hgc_to_hgs(hgccoord, hgsframe):
    """
    Convert from Heliograpic Carrington to Heliographic Stonyhurst.
    """
    if hgsframe.obstime is None or np.any(hgsframe.obstime != hgccoord.obstime):
        raise ValueError("Can not transform from Heliographic Carrington to "
                         "Heliographic Stonyhurst, unless both frames have matching obstime.")
    obstime = hgsframe.obstime
    s_lon = hgccoord.spherical.lon - _carrington_offset(obstime).to(
        u.deg)
    representation = SphericalRepresentation(s_lon, hgccoord.spherical.lat,
                                             hgccoord.spherical.distance)

    return hgsframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Heliocentric,
                                 Helioprojective)
def hcc_to_hpc(helioccoord, heliopframe):
    """
    Convert from Heliocentic Cartesian to Helioprojective Cartesian.
    """
    if not _observers_are_equal(helioccoord.observer, heliopframe.observer):
        heliocframe = Heliocentric(observer=heliopframe.observer)
        new_helioccoord = helioccoord.transform_to(heliocframe)
        helioccoord = new_helioccoord

    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    # d is calculated as the distance between the points
    # (x,y,z) and (0,0,D0).
    distance = np.sqrt(x**2 + y**2 + (helioccoord.observer.radius - z)**2)

    hpcx = np.rad2deg(np.arctan2(x, helioccoord.observer.radius - z))
    hpcy = np.rad2deg(np.arcsin(y / distance))

    representation = SphericalRepresentation(hpcx, hpcy,
                                             distance.to(u.km))

    return heliopframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Helioprojective,
                                 Heliocentric)
def hpc_to_hcc(heliopcoord, heliocframe):
    """
    Convert from Helioprojective Cartesian to Heliocentric Cartesian.
    """
    if not _observers_are_equal(heliopcoord.observer, heliocframe.observer):
        heliocframe_heliopobs = Heliocentric(observer=heliopcoord.observer)
        helioccoord_heliopobs = heliopcoord.transform_to(heliocframe_heliopobs)
        helioccoord = helioccoord_heliopobs.transform_to(heliocframe)
        return helioccoord

    if not isinstance(heliopcoord.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform helioprojective coordinates to "
                           "heliocentric coordinates for observer '{}' "
                           "without `obstime` being specified.".format(heliopcoord.observer))

    heliopcoord = heliopcoord.calculate_distance()
    x = np.deg2rad(heliopcoord.Tx)
    y = np.deg2rad(heliopcoord.Ty)

    cosx = np.cos(x)
    sinx = np.sin(x)
    cosy = np.cos(y)
    siny = np.sin(y)

    rx = (heliopcoord.distance.to(u.m)) * cosy * sinx
    ry = (heliopcoord.distance.to(u.m)) * siny
    rz = (heliopcoord.observer.radius.to(u.m)) - (
        heliopcoord.distance.to(u.m)) * cosy * cosx

    representation = CartesianRepresentation(
        rx.to(u.km), ry.to(u.km), rz.to(u.km))
    return heliocframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Heliocentric,
                                 HeliographicStonyhurst)
def hcc_to_hgs(helioccoord, heliogframe):
    """
    Convert from Heliocentric Cartesian to Heliographic Stonyhurst.
    """
    if not isinstance(helioccoord.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform heliocentric coordinates to "
                           "heliographic coordinates for observer '{}' "
                           "without `obstime` being specified.".format(helioccoord.observer))

    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    l0_rad = helioccoord.observer.lon
    b0_deg = helioccoord.observer.lat

    cosb = np.cos(np.deg2rad(b0_deg))
    sinb = np.sin(np.deg2rad(b0_deg))

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + l0_rad
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)

    representation = SphericalRepresentation(
        np.rad2deg(hgln), np.rad2deg(hglt), hecr.to(u.km))
    return heliogframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst,
                                 Heliocentric)
def hgs_to_hcc(heliogcoord, heliocframe):
    """
    Convert from Heliographic Stonyhurst to Heliocentric Cartesian.
    """
    hglon = heliogcoord.spherical.lon
    hglat = heliogcoord.spherical.lat
    r = heliogcoord.spherical.distance
    if r.unit is u.one and u.allclose(r, 1*u.one):
        r = np.ones_like(r)
        r *= RSUN_METERS

    if not isinstance(heliocframe.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform heliographic coordinates to "
                           "heliocentric coordinates for observer '{}' "
                           "without `obstime` being specified.".format(heliocframe.observer))

    l0_rad = heliocframe.observer.lon.to(u.rad)
    b0_deg = heliocframe.observer.lat

    lon = np.deg2rad(hglon)
    lat = np.deg2rad(hglat)

    cosb = np.cos(b0_deg.to(u.rad))
    sinb = np.sin(b0_deg.to(u.rad))

    lon = lon - l0_rad

    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)

    x = r * cosy * sinx
    y = r * (siny * cosb - cosy * cosx * sinb)
    zz = r * (siny * sinb + cosy * cosx * cosb)

    representation = CartesianRepresentation(
        x.to(u.km), y.to(u.km), zz.to(u.km))

    return heliocframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Helioprojective,
                                 Helioprojective)
def hpc_to_hpc(heliopcoord, heliopframe):
    """
    This converts from HPC to HPC, with different observer location parameters.
    It does this by transforming through HGS.
    """
    if (heliopcoord.observer == heliopframe.observer or
        (u.allclose(heliopcoord.observer.lat, heliopframe.observer.lat) and
         u.allclose(heliopcoord.observer.lon, heliopframe.observer.lon) and
         u.allclose(heliopcoord.observer.radius, heliopframe.observer.radius))):
        return heliopframe.realize_frame(heliopcoord._data)

    if not isinstance(heliopframe.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform between helioprojective frames "
                           "without `obstime` being specified for observer {}.".format(heliopframe.observer))
    if not isinstance(heliopcoord.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform between helioprojective frames "
                           "without `obstime` being specified for observer {}.".format(heliopcoord.observer))

    hgs = heliopcoord.transform_to(HeliographicStonyhurst)
    hgs.observer = heliopframe.observer
    hpc = hgs.transform_to(heliopframe)

    return hpc


def _rotation_matrix_repr_to_repr(start_representation, end_representation):
    """
    Return the matrix for the direct rotation from one representation to a second representation.
    The representations need not be normalized first.
    """
    A = start_representation.to_cartesian()
    B = end_representation.to_cartesian()
    rotation_axis = A.cross(B)
    rotation_angle = -np.arccos(A.dot(B) / (A.norm() * B.norm()))  # negation is required

    # This line works around some input/output quirks of Astropy's rotation_matrix()
    matrix = np.array(rotation_matrix(rotation_angle, rotation_axis.xyz.value.tolist()))
    return matrix


def _rotation_matrix_reprs_to_xz_about_z(representations):
    """
    Return one or more matrices for rotating one or more representations around the Z axis into the
    XZ plane.
    """
    A = representations.to_cartesian()

    # Zero out the Z components
    # (The additional transpose operations are to handle both scalar and array inputs)
    A_no_z = CartesianRepresentation((A.xyz.T * [1, 1, 0]).T)

    # Rotate the resulting vector to the X axis
    x_axis = CartesianRepresentation(1, 0, 0)
    if A_no_z.isscalar:
        matrix = _rotation_matrix_repr_to_repr(A_no_z, x_axis)
    else:
        matrix_list = [_rotation_matrix_repr_to_repr(vect, x_axis) for vect in A_no_z]
        matrix = np.stack(matrix_list)

    return matrix


def _sun_earth_icrf(time):
    """
    Return the Sun-Earth vector for ICRF-based frames.
    """
    sun_pos_icrs = get_body_barycentric('sun', time)
    earth_pos_icrs = get_body_barycentric('earth', time)
    return earth_pos_icrs - sun_pos_icrs


# The Sun's north pole is oriented RA=286.13 deg, dec=63.87 deg in ICRS, and thus HCRS as well
# (See Archinal et al. 2011,
#   "Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009")
# The orientation of the north pole in ICRS/HCRS is assumed to be constant in time
_SOLAR_NORTH_POLE_HCRS = UnitSphericalRepresentation(lon=286.13*u.deg, lat=63.87*u.deg)


# Calculate the rotation matrix to de-tilt the Sun's rotation axis to be parallel to the Z axis
_SUN_DETILT_MATRIX = _rotation_matrix_repr_to_repr(_SOLAR_NORTH_POLE_HCRS,
                                                   CartesianRepresentation(0, 0, 1))


@frame_transform_graph.transform(DynamicMatrixTransform, HCRS, HeliographicStonyhurst)
def hcrs_to_hgs(hcrscoord, hgsframe):
    """
    Convert from HCRS to Heliographic Stonyhurst (HGS).

    HGS shares the same origin (the Sun) as HCRS, but has its Z axis aligned with the Sun's
    rotation axis and its X axis aligned with the projection of the Sun-Earth vector onto the Sun's
    equatorial plane (i.e., the component of the Sun-Earth vector perpendicular to the Z axis).
    Thus, the transformation matrix is the product of the matrix to align the Z axis (by de-tilting
    the Sun's rotation axis) and the matrix to align the X axis.  The first matrix is independent
    of time and is pre-computed, while the second matrix depends on the time-varying Sun-Earth
    vector.
    """
    if hgsframe.obstime is None:
        raise ValueError("To perform this transformation the coordinate"
                         " Frame needs an obstime Attribute")

    # Determine the Sun-Earth vector in HCRS
    sun_earth = _sun_earth_icrf(hgsframe.obstime)

    # De-tilt the Sun-Earth vector to the frame with the Sun's rotation axis parallel to the Z axis
    sun_earth_detilt = sun_earth.transform(_SUN_DETILT_MATRIX)

    # Rotate the Sun-Earth vector about the Z axis so that it lies in the XZ plane
    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(sun_earth_detilt)

    return matrix_product(rot_matrix, _SUN_DETILT_MATRIX)


@frame_transform_graph.transform(DynamicMatrixTransform, HeliographicStonyhurst, HCRS)
def hgs_to_hcrs(hgscoord, hcrsframe):
    """
    Convert from Heliographic Stonyhurst to HCRS.
    """
    return matrix_transpose(hcrs_to_hgs(hcrsframe, hgscoord))


@frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst, HeliographicStonyhurst)
def hgs_to_hgs(from_coo, to_frame):
    """
    Convert between two Heliographic Stonyhurst frames.
    """
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HCRS).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransform, HeliographicCarrington, HeliographicCarrington)
def hgc_to_hgc(from_coo, to_frame):
    """
    Convert between two Heliographic Carrington frames.
    """
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HeliographicStonyhurst(obstime=from_coo.obstime)).\
               transform_to(HeliographicStonyhurst(obstime=to_frame.obstime)).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransform, Heliocentric, Heliocentric)
def hcc_to_hcc(hcccoord, hccframe):
    """
    Convert from  Heliocentric to Heliocentric
    """
    if _observers_are_equal(hcccoord.observer, hccframe.observer, string_ok=True):
        return hccframe.realize_frame(hcccoord._data)

    hgscoord = hcccoord.transform_to(HeliographicStonyhurst)
    hgscoord.observer = hccframe.observer
    hcccoord = hgscoord.transform_to(hccframe)

    return hcccoord


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricMeanEcliptic, HeliocentricEarthEcliptic)
def hme_to_hee(hmecoord, heeframe):
    """
    Convert from Heliocentric Mean Ecliptic to Heliocentric Earth Ecliptic
    """
    # Convert to the HME frame with mean equinox of date at the HEE obstime, through HCRS
    int_frame = HeliocentricMeanEcliptic(obstime=heeframe.obstime, equinox=heeframe.obstime)
    int_coord = hmecoord.transform_to(HCRS).transform_to(int_frame)

    # Get the Sun-Earth vector in the new HME frame
    sun_earth = HCRS(_sun_earth_icrf(int_coord.obstime), obstime=int_coord.obstime)
    sun_earth_hme = sun_earth.transform_to(int_coord).cartesian

    # Rotate the Sun-Earth vector about the Z axis so that it lies in the XZ plane
    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(sun_earth_hme)

    newrepr = int_coord.cartesian.transform(rot_matrix)
    return heeframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricEarthEcliptic, HeliocentricMeanEcliptic)
def hee_to_hme(heecoord, hmeframe):
    """
    Convert from Heliocentric Earth Ecliptic to Heliocentric Mean Ecliptic
    """
    # Get the Sun-Earth vector in the HME frame with mean equinox of date
    int_frame = HeliocentricMeanEcliptic(obstime=heecoord.obstime, equinox=heecoord.obstime)
    sun_earth = HCRS(_sun_earth_icrf(int_frame.obstime), obstime=int_frame.obstime)
    sun_earth_int = sun_earth.transform_to(int_frame).cartesian

    # Rotate the Sun-Earth vector about the Z axis so that it lies in the XZ plane
    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(sun_earth_int)

    # Reverse the rotation to convert from HEE to the HME frame with mean equinox of date
    int_repr = heecoord.cartesian.transform(matrix_transpose(rot_matrix))
    int_coord = int_frame.realize_frame(int_repr)

    return int_coord.transform_to(HCRS).transform_to(hmeframe)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricEarthEcliptic, HeliocentricEarthEcliptic)
def hee_to_hee(from_coo, to_frame):
    """
    Convert between two Heliocentric Earth Ecliptic frames.
    """
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HCRS).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricEarthEcliptic, GeocentricSolarEcliptic)
def hee_to_gse(heecoord, gseframe):
    """
    Convert from Heliocentric Earth Ecliptic to Geocentric Solar Ecliptic
    """
    # Use an intermediate frame of HEE at the GSE observation time
    int_frame = HeliocentricEarthEcliptic(obstime=gseframe.obstime)
    int_coord = heecoord.transform_to(int_frame)

    # Get the Sun-Earth vector in the intermediate frame
    sun_earth = HCRS(_sun_earth_icrf(int_frame.obstime), obstime=int_frame.obstime)
    sun_earth_int = sun_earth.transform_to(int_frame).cartesian

    # Find the Earth-object vector in the intermediate frame
    earth_object_int = int_coord.cartesian - sun_earth_int

    # Flip the vector in X and Y, but leave Z untouched
    # (The additional transpose operations are to handle both scalar and array inputs)
    newrepr = CartesianRepresentation((earth_object_int.xyz.T * [-1, -1, 1]).T)

    return gseframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GeocentricSolarEcliptic, HeliocentricEarthEcliptic)
def gse_to_hee(gsecoord, heeframe):
    """
    Convert from Geocentric Solar Ecliptic to Heliocentric Earth Ecliptic
    """
    # Use an intermediate frame of HEE at the GSE observation time
    int_frame = HeliocentricEarthEcliptic(obstime=gsecoord.obstime)

    # Get the Sun-Earth vector in the intermediate frame
    sun_earth = HCRS(_sun_earth_icrf(int_frame.obstime), obstime=int_frame.obstime)
    sun_earth_int = sun_earth.transform_to(int_frame).cartesian

    # Find the Earth-object vector in the intermediate frame
    # Flip the vector in X and Y, but leave Z untouched
    # (The additional transpose operations are to handle both scalar and array inputs)
    earth_object_int = CartesianRepresentation((gsecoord.cartesian.xyz.T * [-1, -1, 1]).T)

    # Find the Sun-object vector in the intermediate frame
    sun_object_int = sun_earth_int + earth_object_int
    int_coord = int_frame.realize_frame(sun_object_int)

    return int_coord.transform_to(heeframe)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GeocentricSolarEcliptic, GeocentricSolarEcliptic)
def gse_to_gse(from_coo, to_frame):
    """
    Convert between two Geocentric Solar Ecliptic frames.
    """
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HeliocentricEarthEcliptic).transform_to(to_frame)


def _rotation_matrix_hme_to_hci(hmeframe):
    """
    Return the rotation matrix from HME to HCI at the same observation time
    """
    # Get the ecliptic pole and the solar rotation axis
    ecliptic_pole = hmeframe.realize_frame(CartesianRepresentation(0, 0, 1)*u.m)
    solar_rot_axis = HeliographicStonyhurst(CartesianRepresentation(0, 0, 1)*u.m,
                                            obstime=hmeframe.obstime).transform_to(hmeframe)

    # Align the solar rotation axis with the Z axis
    detilt_matrix = _rotation_matrix_repr_to_repr(solar_rot_axis.cartesian,
                                                  CartesianRepresentation(0, 0, 1))
    detilted_ecliptic_pole = ecliptic_pole.cartesian.transform(detilt_matrix)

    # Then align the de-tilted ecliptic pole with the Y axis, which aligns the solar ascending node
    # with the X axis
    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(detilted_ecliptic_pole)
    x_to_y_matrix = _rotation_matrix_repr_to_repr(CartesianRepresentation(1, 0, 0),
                                                  CartesianRepresentation(0, 1, 0))
    rot_matrix = matrix_product(x_to_y_matrix, rot_matrix)

    return matrix_product(rot_matrix, detilt_matrix)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricMeanEcliptic, HeliocentricInertial)
def hme_to_hci(hmecoord, hciframe):
    """
    Convert from Heliocentric Mean Ecliptic to Heliocentric Inertial
    """
    # Convert to the HME frame with mean J2000.0 ecliptic at the HEE obstime, through HCRS
    int_frame = HeliocentricMeanEcliptic(obstime=hciframe.obstime, equinox='J2000.0')
    int_coord = hmecoord.transform_to(HCRS).transform_to(int_frame)

    # Rotate the intermediate coord to the HEE frame
    total_matrix = _rotation_matrix_hme_to_hci(int_frame)
    newrepr = int_coord.cartesian.transform(total_matrix)

    return hciframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricInertial, HeliocentricMeanEcliptic)
def hci_to_hme(hcicoord, hmeframe):
    """
    Convert from Heliocentric Inertial to Heliocentric Mean Ecliptic
    """
    # Use the intermediate frame of HME with mean J2000.0 ecliptic at HCI obstime
    int_frame = HeliocentricMeanEcliptic(obstime=hcicoord.obstime, equinox='J2000.0')

    # Convert the HCI coord to the intermediate frame
    total_matrix = matrix_transpose(_rotation_matrix_hme_to_hci(int_frame))
    newrepr = hcicoord.cartesian.transform(total_matrix)
    int_coord = int_frame.realize_frame(newrepr)

    # Convert to the final frame through HCRS
    return int_coord.transform_to(HCRS).transform_to(hmeframe)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricInertial, HeliocentricInertial)
def hci_to_hci(from_coo, to_frame):
    """
    Convert between two Heliocentric Inertial frames.
    """
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HCRS).transform_to(to_frame)


def _make_sunpy_graph():
    """
    Culls down the full transformation graph for SunPy purposes and returns the string version
    """
    # Frames to keep in the transformation graph
    keep_list = ['icrs', 'hcrs', 'heliocentrictrueecliptic', 'heliocentricmeanecliptic',
                 'heliographic_stonyhurst', 'heliographic_carrington',
                 'heliocentric', 'helioprojective',
                 'heliocentricearthecliptic', 'geocentricsolarecliptic',
                 'heliocentricinertial',
                 'gcrs', 'precessedgeocentric', 'geocentrictrueecliptic', 'geocentricmeanecliptic',
                 'cirs', 'altaz', 'itrs']

    global frame_transform_graph
    backup_graph = deepcopy(frame_transform_graph)

    small_graph = deepcopy(frame_transform_graph)
    cull_list = [name for name in small_graph.get_names() if name not in keep_list]
    cull_frames = [small_graph.lookup_name(name) for name in cull_list]

    for frame in cull_frames:
        # Remove the part of the graph where the unwanted frame is the source frame
        if frame in small_graph._graph:
            del small_graph._graph[frame]

        # Remove all instances of the unwanted frame as the destination frame
        for entry in small_graph._graph:
            if frame in small_graph._graph[entry]:
                del (small_graph._graph[entry])[frame]

    # Clean up the node list
    for name in cull_list:
        small_graph._cached_names.pop(name)

    _add_astropy_node(small_graph)

    # Overwrite the main transform graph
    frame_transform_graph = small_graph

    docstr = make_transform_graph_docs()

    # Restore the main transform graph
    frame_transform_graph = backup_graph

    # Make adjustments to the graph
    docstr = _tweak_graph(docstr)

    return docstr


def _add_astropy_node(graph):
    """
    Add an 'Astropy' node that links to an ICRS node in the graph
    """
    class Astropy(BaseCoordinateFrame):
        name = "REPLACE"

    @graph.transform(FunctionTransform, Astropy, ICRS)
    def fake_transform1():
        pass

    @graph.transform(FunctionTransform, ICRS, Astropy)
    def fake_transform2():
        pass


def _tweak_graph(docstr):
    # Remove Astropy's diagram description
    output = docstr[docstr.find('.. Wrap the graph'):]

    # Change the Astropy node
    output = output.replace('Astropy [shape=oval label="Astropy\\n`REPLACE`"]',
                            'Astropy [shape=box3d style=filled fillcolor=lightcyan '
                            'label="Other frames\\nin Astropy"]')

    # Change the Astropy<->ICRS links to black
    output = output.replace('ICRS -> Astropy[  color = "#783001" ]',
                            'ICRS -> Astropy[  color = "#000000" ]')
    output = output.replace('Astropy -> ICRS[  color = "#783001" ]',
                            'Astropy -> ICRS[  color = "#000000" ]')

    # Set the nodes to be filled and cyan by default
    output = output.replace('AstropyCoordinateTransformGraph {',
                            'AstropyCoordinateTransformGraph {\n'
                            '        node [style=filled fillcolor=lightcyan]')

    # Set the nodes for SunPy frames to be white
    sunpy_frames = ['HeliographicStonyhurst', 'HeliographicCarrington',
                    'Heliocentric', 'Helioprojective',
                    'HeliocentricEarthEcliptic', 'GeocentricSolarEcliptic',
                    'HeliocentricInertial']
    for frame in sunpy_frames:
        output = output.replace(frame + ' [', frame + ' [fillcolor=white ')

    output = output.replace('<ul>\n\n',
                            '<ul>\n\n' +
                            _add_legend_row('SunPy frames', 'white') +
                            _add_legend_row('Astropy frames', 'lightcyan'))

    return output


def _add_legend_row(label, color):
    row = '        <li style="list-style: none;">\n'\
          '            <p style="font-size: 12px;line-height: 24px;font-weight: normal;'\
          'color: #848484;padding: 0;margin: 0;">\n'\
          '                <b>' + label + ':</b>\n'\
          '                    <span class="dot" style="height: 20px;width: 40px;'\
          'background-color: ' + color + ';border-radius: 50%;border: 1px solid black;'\
          'display: inline-block;"></span>\n'\
          '            </p>\n'\
          '        </li>\n\n\n'
    return row
