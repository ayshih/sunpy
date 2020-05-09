"""
Coordinate frames that are defined relative to other frames
"""

from astropy import units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import ConvertError, ICRS, SkyCoord, SphericalRepresentation
from astropy.coordinates.attributes import Attribute, QuantityAttribute
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransform
import numpy as np

from sunpy import log
from sunpy.time import parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring

from .frameattributes import TimeFrameAttributeSunPy
from .frames import _J2000, HeliocentricInertial, SunPyBaseCoordinateFrame
from .offset_frame import NorthOffsetFrame
from .transformations import _transformation_debug

__all__ = ['NorthOffsetFrame', 'RotatedSunFrame', 'AstrometricFrame']


# The code for NorthOffsetFrame currently lives in `offset_frame.py`
# This changes its module to this file so that the docs think that the class is local
NorthOffsetFrame.__module__ = __name__


_rotatedsun_cache = {}
_astrometric_cache = {}

# The threshold distance that distinguishes between solar-system bodies and cosmic objects
_DISTANCE_THRESHOLD = 1*u.lyr


def _conform_frame_input(frame):
    # If a SkyCoord is provided, use the underlying frame
    if hasattr(frame, 'frame'):
        frame = frame.frame

    # The frame needs to be a SunPy frame to have the overridden size() property
    # This can be removed with Astropy 4.0
    if not isinstance(frame, SunPyBaseCoordinateFrame):
        raise TypeError("Only SunPy coordinate frames are currently supported.")

    return frame


def _make_rotatedsun_cls(framecls):
    """
    Create a new class that is the rotated-Sun frame for a specific class of
    base frame. If such a class has already been created for this frame, the
    same class will be returned.

    This function is necessary because frame transformations depend
    on connection between specific frame *classes*.  So each type of frame
    needs its own distinct rotated-Sun frame class.  This function generates
    just that class, as well as ensuring that only one example of such a class
    actually gets created in any given Python session.
    """
    # This code reuses significant code from Astropy's implementation of SkyOffsetFrame
    # See licenses/ASTROPY.rst

    if framecls in _rotatedsun_cache:
        return _rotatedsun_cache[framecls]

    # Obtain the base frame's metaclass by getting the type of the base frame's class
    framemeta = type(framecls)

    # Subclass the metaclass for the RotatedSunFrame from the base frame's metaclass
    class RotatedSunMeta(framemeta):
        """
        This metaclass renames the class to be "RotatedSun<framecls>".
        """
        def __new__(cls, name, bases, members):
            newname = name[:-5] if name.endswith('Frame') else name
            newname += framecls.__name__

            return super().__new__(cls, newname, bases, members)

    # We need this to handle the intermediate metaclass correctly, otherwise we could
    # just subclass RotatedSunFrame.
    _RotatedSunFramecls = RotatedSunMeta('RotatedSunFrame',
                                         (RotatedSunFrame, framecls),
                                         {'__doc__': RotatedSunFrame.__doc__})

    @frame_transform_graph.transform(FunctionTransform, _RotatedSunFramecls, _RotatedSunFramecls)
    @_transformation_debug(f"{_RotatedSunFramecls.__name__}->{_RotatedSunFramecls.__name__}")
    def rotatedsun_to_rotatedsun(from_rotatedsun_coord, to_rotatedsun_frame):
        """Transform between two rotated-Sun frames."""
        # This transform goes through the parent frames on each side.
        # from_frame -> from_frame.base -> to_frame.base -> to_frame
        intermediate_from = from_rotatedsun_coord.transform_to(from_rotatedsun_coord.base)
        intermediate_to = intermediate_from.transform_to(to_rotatedsun_frame.base)
        return intermediate_to.transform_to(to_rotatedsun_frame)

    @frame_transform_graph.transform(FunctionTransform, framecls, _RotatedSunFramecls)
    @_transformation_debug(f"{framecls.__name__}->{_RotatedSunFramecls.__name__}")
    def reference_to_rotatedsun(reference_coord, rotatedsun_frame):
        # Transform to HCI
        hci_frame = HeliocentricInertial(obstime=rotatedsun_frame.base.obstime)
        hci_coord = reference_coord.transform_to(hci_frame)
        oldrepr = hci_coord.spherical

        # Rotate the coordinate in HCI
        from sunpy.physics.differential_rotation import diff_rot
        log.debug(f"Applying {rotatedsun_frame.duration} of solar rotation")
        newlon = oldrepr.lon - diff_rot(rotatedsun_frame.duration,
                                        oldrepr.lat,
                                        rot_type=rotatedsun_frame.rotation_model,
                                        frame_time='sidereal')
        newrepr = SphericalRepresentation(newlon, oldrepr.lat, oldrepr.distance)

        # Transform back from HCI
        new_coord = hci_coord.realize_frame(newrepr).transform_to(rotatedsun_frame.base)
        return rotatedsun_frame.realize_frame(new_coord.data)

    @frame_transform_graph.transform(FunctionTransform, _RotatedSunFramecls, framecls)
    @_transformation_debug(f"{_RotatedSunFramecls.__name__}->{framecls.__name__}")
    def rotatedsun_to_reference(rotatedsun_coord, reference_frame):
        # Transform to HCI
        from_coord = rotatedsun_coord.base.realize_frame(rotatedsun_coord.data)
        hci_coord = from_coord.transform_to(HeliocentricInertial(obstime=reference_frame.obstime))
        oldrepr = hci_coord.spherical

        # Rotate the coordinate in HCI
        from sunpy.physics.differential_rotation import diff_rot
        log.debug(f"Applying {rotatedsun_coord.duration} of solar rotation")
        newlon = oldrepr.lon + diff_rot(rotatedsun_coord.duration,
                                        oldrepr.lat,
                                        rot_type=rotatedsun_coord.rotation_model,
                                        frame_time='sidereal')
        newrepr = SphericalRepresentation(newlon, oldrepr.lat, oldrepr.distance)

        # Transform back from HCI
        hci_coord = HeliocentricInertial(newrepr, obstime=reference_frame.obstime)
        return hci_coord.transform_to(reference_frame)

    _rotatedsun_cache[framecls] = _RotatedSunFramecls
    return _RotatedSunFramecls


@add_common_docstring(**_variables_for_parse_time_docstring())
class RotatedSunFrame:
    """
    A frame that applies solar rotation to a base coordinate frame.

    .. note::

        See :ref:`sunpy-coordinates-rotatedsunframe` for how to use this class.

    In essence, the coordinate axes of the frame are distorted by differential solar rotation.
    This allows using a coordinate representation at one time (at the ``obstime`` of the base
    coordinate frame) to point to a location at a different time that has been differentially
    rotated by the time difference (``duration``).

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or ``None``
        A representation object or ``None`` to have no data.  Alternatively, use coordinate
        component keyword arguments, which depend on the base frame.
    base : `~astropy.coordinates.SkyCoord` or low-level coordinate object.
        The coordinate which specifies the base coordinate frame.  The frame must be a SunPy frame.
    duration : `~astropy.units.Quantity`
        The duration of solar rotation (defaults to zero days).
    rotated_time : {parse_time_types}
        The time to rotate the Sun to.  If provided, ``duration`` will be set to the difference
        between this time and the observation time in ``base``.
    rotation_model : `str`
        Accepted model names are ``'howard'`` (default), ``'snodgrass'``, and ``'allen'``.

    Notes
    -----
    ``RotatedSunFrame`` is a factory class.  That is, the objects that it
    yields are *not* actually objects of class ``RotatedSunFrame``.  Instead,
    distinct classes are created on-the-fly for whatever the frame class is
    of ``base``.
    """
    # This code reuses significant code from Astropy's implementation of SkyOffsetFrame
    # See licenses/ASTROPY.rst

    # Even though the frame attribute `base` is a coordinate frame, we use `Attribute` instead of
    # `CoordinateAttribute` because we are preserving the supplied frame rather than converting to
    # a common frame.
    base = Attribute()

    duration = QuantityAttribute(default=0*u.day)
    rotation_model = Attribute(default='howard')

    def __new__(cls, *args, **kwargs):
        # We don't want to call this method if we've already set up
        # an rotated-Sun frame for this class.
        if not (issubclass(cls, RotatedSunFrame) and cls is not RotatedSunFrame):
            # We get the base argument, and handle it here.
            base_frame = kwargs.get('base', None)
            if base_frame is None:
                raise TypeError("Can't initialize a RotatedSunFrame without a `base` keyword.")
            kwargs['base'] = _conform_frame_input(base_frame)

            newcls = _make_rotatedsun_cls(kwargs['base'].__class__)
            return newcls.__new__(newcls, *args, **kwargs)

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)
        return super().__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        # Validate inputs
        if kwargs['base'].obstime is None:
            raise ValueError("The base coordinate frame must have a defined `obstime`.")

        if 'rotated_time' in kwargs:
            rotated_time = parse_time(kwargs['rotated_time'])
            kwargs['duration'] = (rotated_time - kwargs['base'].obstime).to('day')
            kwargs.pop('rotated_time')

        super().__init__(*args, **kwargs)

        # Move data out from the base frame
        if self.base.has_data:
            if not self.has_data:
                # If the duration is an array but the data is scalar, upgrade data to an array
                if self.base.data.isscalar and not self.duration.isscalar:
                    self._data = self.base.data._apply('repeat', self.duration.shape)
                else:
                    self._data = self.base.data
            self._base = self.base.replicate_without_data()

    def as_base(self):
        """
        Returns a coordinate with the current representation and in the base coordinate frame.

        This method can be thought of as "removing" the
        `~sunpy.coordinates.metaframes.RotatedSunFrame` layer.  Be aware that this method is not
        merely a coordinate transformation, because this method changes the location in inertial
        space that is being pointed to.
        """
        return self.base.realize_frame(self.data)

    @property
    def rotated_time(self):
        """
        Returns the sum of the base frame's observation time and the rotation of duration.
        """
        return self.base.obstime + self.duration


def _make_astrometric_cls(framecls):
    """
    Create a new class that is the astrometric frame for a specific class of
    base frame. If such a class has already been created for this frame, the
    same class will be returned.

    This function is necessary because frame transformations depend
    on connection between specific frame *classes*.  So each type of frame
    needs its own distinct astrometric frame class.  This function generates
    just that class, as well as ensuring that only one example of such a class
    actually gets created in any given Python session.
    """
    # This code reuses significant code from Astropy's implementation of SkyOffsetFrame
    # See licenses/ASTROPY.rst

    if framecls in _astrometric_cache:
        return _astrometric_cache[framecls]

    # Obtain the base frame's metaclass by getting the type of the base frame's class
    framemeta = type(framecls)

    # Subclass the metaclass for the AstrometricFrame from the base frame's metaclass
    class AstrometricMeta(framemeta):
        """
        This metaclass renames the class to be "Astrometric<framecls>".
        """
        def __new__(cls, name, bases, members):
            newname = name[:-5] if name.endswith('Frame') else name
            newname += framecls.__name__

            return super().__new__(cls, newname, bases, members)

    # We need this to handle the intermediate metaclass correctly, otherwise we could
    # just subclass AstrometricFrame.
    _AstrometricFramecls = AstrometricMeta('AstrometricFrame',
                                         (AstrometricFrame, framecls),
                                         {'__doc__': AstrometricFrame.__doc__})

    @frame_transform_graph.transform(FunctionTransform, _AstrometricFramecls, _AstrometricFramecls)
    @_transformation_debug(f"{_AstrometricFramecls.__name__}->{_AstrometricFramecls.__name__}")
    def astrometric_to_astrometric(from_astrometric_coord, to_astrometric_frame):
        """Transform between two astrometric frames."""
        # This transform goes through the parent frames on each side.
        # from_frame -> from_frame.base -> to_frame.base -> to_frame
        intermediate_from = from_astrometric_coord.transform_to(from_astrometric_coord.base)
        intermediate_to = intermediate_from.transform_to(to_astrometric_frame.base)
        return intermediate_to.transform_to(to_astrometric_frame)

    @frame_transform_graph.transform(FunctionTransform, framecls, _AstrometricFramecls)
    @_transformation_debug(f"{framecls.__name__}->{_AstrometricFramecls.__name__}")
    def reference_to_astrometric(reference_coord, astrometric_frame):
        reference_coord = reference_coord.transform_to(astrometric_frame.base)
        icrs_coord = reference_coord.transform_to(ICRS)

        posvel = icrs_coord.cartesian
        if 's' in posvel.differentials:
            distance = posvel.norm()
            if np.all(distance < _DISTANCE_THRESHOLD):  # solar-system body
                pos = posvel.without_differentials()
                vel = posvel.differentials['s'].to_cartesian()

                # The light travel time (dt) for a body in linear motion satisfies the equation:
                #   norm(D - V * dt) = c * dt
                # where D is the observer-body vector, V is the body velocity vector,
                # and c is the speed of light.  This turns into a quadratic equation in dt, and
                # only the positive solution is valid:
                #   dt = (sqrt((D dot V)^2 + (D dot D) * c2_V2) - D dot V) / c2_V2
                # where c2_V2 = c^2 - V dot V
                D = pos - astrometric_frame._observer_icrs_cartesian
                c2_V2 = speed_of_light ** 2 - vel.dot(vel)
                DdotV = D.dot(vel)
                dt = (np.sqrt(DdotV ** 2 + D.dot(D) * c2_V2) - DdotV) / c2_V2
                log.debug(f"Retarding for {dt.to(u.s)} of light travel time")

                pos -= vel * dt
                posvel = pos.with_differentials(posvel.differentials)
                icrs_coord = ICRS(posvel)
            elif np.all(distance >= _DISTANCE_THRESHOLD):  # cosmic object
                dt = (astrometric_frame.obstime - astrometric_frame.ref_epoch).to(u.yr)
                log.debug(f"Propagating forward by {dt} after {astrometric_frame.ref_epoch}")
                icrs_coord = SkyCoord(icrs_coord).apply_space_motion(dt=dt).frame
            else:  # a mix
                raise ConvertError("The transformation cannot handle a mix of solar-system bodies "
                                   "and cosmic objects.")

        reference_coord = icrs_coord.transform_to(reference_coord)
        return astrometric_frame.realize_frame(reference_coord.data)

    @frame_transform_graph.transform(FunctionTransform, _AstrometricFramecls, framecls)
    @_transformation_debug(f"{_AstrometricFramecls.__name__}->{framecls.__name__}")
    def astrometric_to_reference(astrometric_coord, reference_frame):
        base_coord = astrometric_coord.as_base()
        icrs_coord = base_coord.transform_to(ICRS)

        posvel = icrs_coord.cartesian
        if 's' in posvel.differentials:
            distance = posvel.norm()
            if np.all(distance < _DISTANCE_THRESHOLD):  # solar-system body
                pos = posvel.without_differentials()

                # The light travel time is straightforward to compute for an astrometric location
                dt = (pos - astrometric_coord._observer_icrs_cartesian).norm() / speed_of_light
                log.debug(f"Advancing for {dt.to(u.s)} of light travel time")

                pos += posvel.differentials['s'] * dt
                posvel = pos.with_differentials(posvel.differentials)
                icrs_coord = ICRS(posvel)
            elif np.all(distance >= _DISTANCE_THRESHOLD):  # cosmic object
                dt = (astrometric_coord.obstime - astrometric_coord.ref_epoch).to(u.yr)
                log.debug(f"Propagating backward by {dt} before {astrometric_coord.ref_epoch}")
                icrs_coord = SkyCoord(icrs_coord).apply_space_motion(dt=-dt).frame
            else:  # a mix
                raise ConvertError("The transformation cannot handle a mix of solar-system bodies "
                                   "and cosmic objects.")

        base_coord = icrs_coord.transform_to(base_coord)
        return base_coord.transform_to(reference_frame)

    _astrometric_cache[framecls] = _AstrometricFramecls
    return _AstrometricFramecls


@add_common_docstring(**_variables_for_parse_time_docstring())
class AstrometricFrame:
    """
    An astrometric version of a base coordinate frame.

    .. note::

        See :ref:`sunpy-coordinates-astrometricframe` for how to use this class.

    The astrometric location of a body is the location of the body as measured by an observer,
    which depends on the motion of that body:

    * For a body in the solar system, its coordinate normally represents the instantaneous
      location of that body.  However, due to the finite speed of light, the light that reaches an
      observer at a certain observation time was emitted at an earlier time.  If the body is
      moving, then the observer will measure the body to be where it was at that earlier time.
      This effect is also known as "planetary aberration".

    * For a cosmic object, its coordinate normally represents its catalog location at a reference
      epoch (e.g., J2000.0).  If the body is moving (i.e., has non-zero "proper motion" or "radial
      velocity"), its measured location will evolve over time.  In this case, light travel time is
      already included in the catalog location.

    .. warning::
        The astrometric calculation assumes that the body is moving in a straight line with its
        specified velocity.  The accuracy of this assumption depends on the body's true motion and
        the observer-body distance (or light travel time).

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or ``None``
        A representation object or ``None`` to have no data.  Alternatively, use coordinate
        component keyword arguments, which depend on the base frame.
    base : `~astropy.coordinates.SkyCoord` or low-level coordinate object.
        The coordinate which specifies the base coordinate frame.  The frame must be a SunPy frame.
    observer : `~astropy.coordinates.SkyCoord` or low-level coordinate object.
        The coordinate which specifies the observer coordinate frame.  If an observer is not
        provided this way, the observer from ``base`` will be used if it exists.  If an observer is
        not provided in either fashion, astrometric calculations for solar-system bodies cannot be
        performed.  The frame must be a SunPy frame.
    ref_epoch : {parse_time_types}
        The reference epoch for the location of cosmic objects.  Defaults to J2000.0.

    Notes
    -----
    If the coordinate has no velocity included, then its astrometric location is the same as its
    instantaneous location.

    The astrometric calculation distinguishes between bodies in the solar system and cosmic objects
    by the distance from the solar-system barycenter (SSB).  The coordinates of bodies less than 1
    light-year from the SSB are treated as instantaneous locations.  The coordinates of bodies
    greater than 1 light-year from the SSB are treated as locations at the reference epoch.

    As these are "astrometric" coordinates rather than "apparent" coordinates, effects such as
    stellar aberration and gravitational deflection are not included.

    ``AstrometricFrame`` is a factory class.  That is, the objects that it
    yields are *not* actually objects of class ``AstrometricFrame``.  Instead,
    distinct classes are created on-the-fly for whatever the frame class is
    of ``base``.
    """
    # This code reuses significant code from Astropy's implementation of SkyOffsetFrame
    # See licenses/ASTROPY.rst

    # Even though the frame attribute `base` is a coordinate frame, we use `Attribute` instead of
    # `CoordinateAttribute` because we are preserving the supplied frame rather than converting to
    # a common frame.
    base = Attribute()
    observer = Attribute()

    ref_epoch = TimeFrameAttributeSunPy(default=_J2000)

    def __new__(cls, *args, **kwargs):
        # We don't want to call this method if we've already set up
        # an astrometric frame for this class.
        if not (issubclass(cls, AstrometricFrame) and cls is not AstrometricFrame):
            # We get the base argument, and handle it here.
            base_frame = kwargs.get('base', None)
            if base_frame is None:
                raise TypeError("Can't initialize an AstrometricFrame without a `base` keyword.")
            kwargs['base'] = _conform_frame_input(base_frame)

            newcls = _make_astrometric_cls(kwargs['base'].__class__)
            return newcls.__new__(newcls, *args, **kwargs)

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)
        return super().__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        # Validate base obstime
        if kwargs['base'].obstime is None:
            raise ValueError("The base coordinate frame must have a defined `obstime`.")

        # If `observer` isn't provided explicitly, try to pull it from `base`
        if kwargs.get('observer', None) is None:
            kwargs['observer'] = getattr(kwargs['base'], 'observer', None)

        # Validate observer obstime if there is an observer
        if kwargs['observer'] is not None:
            kwargs['observer'] = _conform_frame_input(kwargs['observer'])

            if kwargs['observer'].obstime is None:
                raise ValueError("The observer coordinate frame must have a defined `obstime`.")

            if np.any(kwargs['base'].obstime != kwargs['observer'].obstime):
                raise ValueError("The `obstime` for the base coordinate frame must match the "
                                 "`obstime` for the observer coordinate frame.")

        super().__init__(*args, **kwargs)

        # Pre-compute the Cartesian representation of the observer location converted to ICRS
        if self.observer is not None:
            self._observer_icrs_cartesian = self.observer.transform_to(ICRS).cartesian

        # Validate main obstime
        if self.obstime is not None:
            if np.any(self.obstime != self.base.obstime):
                raise ValueError("The `obstime` for the main coordinate frame must match the "
                                 "`obstime` for the base coordinate frame.")
        else:
            self._obstime = self.base.obstime

        # Move data out from the base frame
        if self.base.has_data:
            self._data = self.base.data
            self._base = self.base.replicate_without_data()

        # Check that any frame data is not dimensionless
        if self.has_data:
            if self._data.norm().unit is u.one:
                raise ValueError("The frame data must have distance units.")

    def as_base(self):
        """
        Returns a coordinate with the current representation and in the base coordinate frame.

        This method can be thought of as "removing" the
        `~sunpy.coordinates.metaframes.AstrometricFrame` layer.  Be aware that this method is not
        merely a coordinate transformation, because this method changes the location in inertial
        space that is being pointed to.
        """
        return self.base.realize_frame(self.data)
