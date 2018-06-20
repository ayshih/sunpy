# -*- coding: utf-8 -*-
"""
Allow representation of queries as logic expressions. This module makes
sure that attributes that are combined using the two logic operations AND (&)
and OR (|) always are in disjunctive normal form, that is, there are only two
levels ­- the first being disjunction and the second being conjunction. In other
words, every combinations of attributes looks like this:
(a AND b AND c) OR (d AND e).

Walkers are used to traverse the tree that results from combining attributes.
They are implemented using sunpy.util.multimethod. Multimethods are functions
that are not assigned to classes but still dispatch by type of one or more
of their arguments. For more information about multimethods, refer to
sunpy.util.multimethod.

Please note that & is evaluated first, so A & B | C is equivalent to
(A & B) | C.
"""
import keyword
from collections import defaultdict, namedtuple

from astropy.utils.misc import isiterable

from sunpy.util.multimethod import MultiMethod
from sunpy.extern.six import iteritems

_ATTR_TUPLE = namedtuple("attr", "name name_long des")


def make_tuple():
    return _ATTR_TUPLE([], [], [])


class AttrMeta(type):
    """
    We want to enable automatic discovery via tab completion of subclasses of Attrs.
    To do this we have to create a metaclass that redefines the methods that Python uses to normally do this.
    But also to enable the registration of attributes on import.
    Which would allow `a.Instrument` to be able to tab complete to `a.Instrument.AIA` or `a.Instrument.HMI` which does not happen without this.
    """

    # The aim is to register Attrs as a namedtuple of lists
    _attr_registry = defaultdict(make_tuple)

    def __getattr__(self, item):
        """
        """
        if item in self._attr_registry[self].name:
            import pdb; pdb.set_trace()
            return self(item)
        else:
            return None

    def __dir__(self):
        """
        """
        return super().__dir__() + self._attr_registry[self].name

#    def __repr__(self):
#        return str()


class Attr(metaclass=AttrMeta):
    """ This is the base for all attributes. """
    def __and__(self, other):
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        if self.collides(other):
            return NotImplemented
        if isinstance(other, AttrAnd):
            return AttrAnd([self] + list(other.attrs))
        return AttrAnd([self, other])

    def __hash__(self):
        return hash(frozenset(iteritems(vars(self))))

    def __or__(self, other):
        # Optimization.
        if self == other:
            return self
        if isinstance(other, AttrOr):
            return AttrOr([self] + list(other.attrs))
        return AttrOr([self, other])

    def collides(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        return dict(vars(self)) == dict(vars(other))

    @classmethod
    def update_values(cls, adict):
        """
        This is how clients will register their `Attrs` with the system.
        The input should be a dictionary.
        Each key should be a type of the attrs you want.
        The value should be a list of pairs.
        The pair should be a tuple of (Name, Description).
        If you do not want to add a description, you can put `"N/A"`.
        We sanitize the name and it becomes an attribute on the attrs class.
        # Do we want to have them use type?

        We have an example here for an Instrument Class.

        Example
        -------
        >>> from sunpy.net import attr, attrs
        >>> attr.Attr.update_values({type(attrs.Instrument()): [('AIA', 'AIA is in Space.'), ('HMI', 'HMI is next to AIA.')]})

        """
        for k, v in adict.items():
            if isiterable(v) and not isinstance(v, str):
                for pair in v:
                    if len(pair) != 2:
                        raise Warning
                    else:
                        # Sanity name, we remove all special characters and make it all lower case
                        name = ''.join(char for char in pair[0] if char.isalnum()).lower()
                        if keyword.iskeyword(name) or not name.isidentifier():
                            raise Warning
                        else:
                            cls._attr_registry[k][0].append(name)
                            cls._attr_registry[k][1].append(pair[0])
                            cls._attr_registry[k][2].append(pair[1])
            else:
                raise NotImplementedError


class DummyAttr(Attr):
    """ Empty attribute. Useful for building up queries. Returns other
    attribute when ORed or ANDed. It can be considered an empty query
    that you can use as an initial value if you want to build up your
    query in a loop.

    So, if we wanted an attr matching all the time intervals between the times
    stored as (from, to) tuples in a list, we could do.

    attr = DummyAttr()
    for from\_, to in times:
        attr |= Time(from\_, to)
    """
    def __and__(self, other):
        return other

    def __or__(self, other):
        return other

    def collides(self, other):
        return False

    def __hash__(self):
        return hash(None)

    def __eq__(self, other):
        return isinstance(other, DummyAttr)


class SimpleAttr(Attr):
    """
    An attribute that only has a single value.

    This type of attribute is not a composite and has a single value such as
    ``Instrument('EIT')``.

    Parameters
    ----------
    value : `object`
       The value for the attribute to hold.
    """
    def __init__(self, value):
        Attr.__init__(self)
        self.value = value

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __repr__(self):
        return "<{cname!s}({val!r})>".format(
            cname=self.__class__.__name__, val=self.value)


class AttrAnd(Attr):
    """ Attribute representing attributes ANDed together. """
    def __init__(self, attrs):
        Attr.__init__(self)
        self.attrs = attrs

    def __and__(self, other):
        if any(other.collides(elem) for elem in self.attrs):
            return NotImplemented
        if isinstance(other, AttrAnd):
            return AttrAnd(self.attrs + other.attrs)
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        return AttrAnd(self.attrs + [other])

    __rand__ = __and__

    def __repr__(self):
        return "<AttrAnd({att!r})>".format(att=self.attrs)

    def __eq__(self, other):
        if not isinstance(other, AttrAnd):
            return False
        return set(self.attrs) == set(other.attrs)

    def __hash__(self):
        return hash(frozenset(self.attrs))

    def collides(self, other):
        return any(elem.collides(other) for elem in self.attrs)


class AttrOr(Attr):
    """ Attribute representing attributes ORed together. """
    def __init__(self, attrs):
        Attr.__init__(self)
        self.attrs = attrs

    def __or__(self, other):
        if isinstance(other, AttrOr):
            return AttrOr(self.attrs + other.attrs)
        return AttrOr(self.attrs + [other])

    __ror__ = __or__

    def __and__(self, other):
        return AttrOr([elem & other for elem in self.attrs])

    __rand__ = __and__

    def __xor__(self, other):
        new = AttrOr([])
        for elem in self.attrs:
            try:
                new |= elem ^ other
            except TypeError:
                pass
        return new

    def __contains__(self, other):
        for elem in self.attrs:
            try:
                if other in elem:
                    return True
            except TypeError:
                pass
        return False

    def __repr__(self):
        return "<AttrOr({att!r})>".format(att=self.attrs)

    def __eq__(self, other):
        if not isinstance(other, AttrOr):
            return False
        return set(self.attrs) == set(other.attrs)

    def __hash__(self):
        return hash(frozenset(self.attrs))

    def collides(self, other):
        return all(elem.collides(other) for elem in self.attrs)


class ValueAttr(Attr):
    def __init__(self, attrs):
        Attr.__init__(self)
        self.attrs = attrs

    def __repr__(self):
        return "<ValueAttr({att!r})>".format(att=self.attrs)

    def __hash__(self):
        return hash(frozenset(self.attrs.iteritems.items()))

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.attrs == other.attrs

    def collides(self, other):
        if not isinstance(other, self.__class__):
            return False
        return any(k in other.attrs for k in self.attrs)


class AttrWalker(object):
    def __init__(self):
        self.applymm = MultiMethod(lambda *a, **kw: (a[1], ))
        self.createmm = MultiMethod(lambda *a, **kw: (a[1], ))

    def add_creator(self, *types):
        def _dec(fun):
            for type_ in types:
                self.createmm.add(fun, (type_, ))
            return fun
        return _dec

    def add_applier(self, *types):
        def _dec(fun):
            for type_ in types:
                self.applymm.add(fun, (type_, ))
            return fun
        return _dec

    def add_converter(self, *types):
        def _dec(fun):
            for type_ in types:
                self.applymm.add(self.cv_apply(fun), (type_, ))
                self.createmm.add(self.cv_create(fun), (type_, ))
            return fun
        return _dec

    def cv_apply(self, fun):
        def _fun(*args, **kwargs):
            args = list(args)
            args[1] = fun(args[1])
            return self.applymm(*args, **kwargs)
        return _fun

    def cv_create(self, fun):
        def _fun(*args, **kwargs):
            args = list(args)
            args[1] = fun(args[1])
            return self.createmm(*args, **kwargs)
        return _fun

    def create(self, *args, **kwargs):
        return self.createmm(self, *args, **kwargs)

    def apply(self, *args, **kwargs):
        return self.applymm(self, *args, **kwargs)

    def super_create(self, *args, **kwargs):
        return self.createmm.super(self, *args, **kwargs)

    def super_apply(self, *args, **kwargs):
        return self.applymm.super(self, *args, **kwargs)


def and_(*args):
    """ Trick operator precedence.

    and_(foo < bar, bar < baz)
    """
    value = DummyAttr()
    for elem in args:
        value &= elem
    return value


def or_(*args):
    """ Trick operator precedence.

    or_(foo < bar, bar < baz)
    """
    value = DummyAttr()
    for elem in args:
        value |= elem
    return value
