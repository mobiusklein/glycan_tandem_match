from glypy import Composition
from brainpy import mass_charge_ratio

from ms_deisotope.utils import dict_proxy
from collections import defaultdict


def _formulastr(composition):
    parts = [
        "%s" % k if v == 1
        else "%s%d" % (k, v)
        for k, v in sorted(composition.items())
    ]
    return ''.join(parts)


@dict_proxy("composition")
class Formula(object):
    """Represent a hashable composition wrapper and a backing object supplying
    other attributes and properties. Used to make handling of objects which posses
    composition-like behavior easier.

    This class implements the Mapping interface, supplying its composition as key-value data.

    Attributes
    ----------
    composition : Composition
        The elemental composition to proxy
    data : object
        An arbitrary object to proxy attribute look ups to.
        Usually the source of :attr:`composition`
    name : str
        The formula for the given composition. Used for hashing and equality testing
    """
    def __init__(self, composition, data):
        self.composition = composition
        self.name = _formulastr(composition)
        self.data = data

    def __getattr__(self, attr):
        return getattr(self.data, attr)

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)


@dict_proxy("composition")
class StructureRecord(object):
    """Represent a single glycan structure and all of the fragments of interest
    for database search

    This class implements the Mapping interface, supplying its composition as key-value data.

    This class implements the `root` and `tree` protocols from `glypy` making it compatible with tree
    traversing functions like `plot` and `subtree_search`

    Attributes
    ----------
    fragment_map : dict
        A dictionary mapping fragment name to :class:`glypy.structure.fragment.Fragment`
        objects derived from :attr:`structure`
    structure : glypy.Glycan
        The structure represented
    composition : Composition
        The total composition of the structure represented.
    """
    def __init__(self, structure, fragment_map=None, max_cleavages=2):
        self.structure = structure
        if fragment_map is None:
            fragment_map = {f.name: f for f in structure.fragments(
                "ABCXYZ", max_cleavages=max_cleavages, traversal_method='index')}
        self.fragment_map = fragment_map

    def filter_fragment_type(self, exclude):
        """Remove all fragments matching the given type in `exclude` from
        :attr:`fragment_map`

        Parameters
        ----------
        exclude : str
            The type of fragment to remove. Should be one of A, B, C, X, Y, or Z
        """
        self.fragment_map = {
            f.name: f for f in self.fragment_map.values()
            if not set(exclude) & set(f.kind)
        }

    def __root__(self):
        return self.structure.root

    def __tree__(self):
        return self.structure

    def __repr__(self):
        return "StructureRecord(%s)" % self.structure.name

    def mass(self, charge=0, charge_carrier=Composition("Na").mass):
        """Calculate the mass of this structure from its elemental composition.

        Parameters
        ----------
        charge : int, optional
            Instead of calculating neutral mass, calculate the mass to charge ratio
            with charge = `charge`
        charge_carrier : float, optional
            The charge carrier to use when charge is non-zero. By default this is a sodium.

        Returns
        -------
        float
        """
        if charge == 0:
            return self.structure.mass()
        return mass_charge_ratio(self.structure.mass(), z=charge, charge_carrier=charge_carrier)

    @property
    def composition(self):
        return self.structure.total_composition()


@dict_proxy("composition")
class CommonComposition(object):
    """A type to represent multiple objects which share the same elemental composition,
    holding a list of all the objects which it shares formulae with.

    This class implements the Mapping interface, supplying its composition as key-value data.

    Attributes
    ----------
    base_data : list
        All the objects which are bundled together.
    composition : Composition
        The composition being represented.
    """
    def __init__(self, composition, base_data):
        self.composition = composition
        self.base_data = base_data

    @property
    def name(self):
        return '|'.join(x.name for x in self.base_data)

    @classmethod
    def aggregate(cls, iterable):
        """Combine an iterable of objects which are hashable and function as
        composition mappings into a list of :class:`CommonComposition` instances such
        that each composition is represented exactly once.

        Parameters
        ----------
        iterable : Iterable
            An iterable of any type of object which is both hashable and provides a Composition-like
            interface

        Returns
        -------
        list
        """
        agg = defaultdict(list)
        for item in iterable:
            agg[item].append(item.data)
        results = []
        for key, collected in agg.items():
            results.append(cls(key, collected))
        return results

    def __repr__(self):
        return "CommonComposition(%s)" % dict(self.composition)


class StructureDatabase(object):
    """A quick-to-search database of :class:`StructureRecord` instances
    stored in memory.

    Implements the Sequence interface, with `__iter__`, `__len__`, and `__getitem__`.

    Attributes
    ----------
    structures : list
        A list of :class:`StructureRecord` instances, sorted by mass
    """
    def __init__(self, structures):
        self.structures = list(structures)
        self.structures.sort(key=lambda x: x.mass())

    def __len__(self):
        return len(self.structures)

    def __iter__(self):
        return iter(self.structures)

    def __getitem__(self, index):
        return self.structures[index]

    def search_binary(self, mass, error_tolerance=1e-6):
        """Search within :attr:`structures` for the index of a structure
        with a mass nearest to `mass`, within `error_tolerance`

        Parameters
        ----------
        mass : float
            The neutral mass to search for
        error_tolerance : float, optional
            The approximate error tolerance to accept

        Returns
        -------
        int
            The index of the structure with the nearest mass
        """
        lo = 0
        hi = len(self)

        while hi != lo:
            mid = (hi + lo) / 2
            x = self[mid]
            err = x.mass() - mass
            if abs(err) <= error_tolerance:
                return mid
            elif (hi - lo) == 1:
                return mid
            elif err > 0:
                hi = mid
            elif err < 0:
                lo = mid

    def search_mass(self, mass, error_tolerance=0.1):
        """Search for the set of all items in :attr:`structures` within `error_tolerance` Da
        of the queried `mass`.

        Parameters
        ----------
        mass : float
            The neutral mass to search for
        error_tolerance : float, optional
            The range of mass errors (in Daltons) to allow

        Returns
        -------
        list
            The list of :class:`StructureRecord` instances which meet the criterion
        """
        lo_mass = mass - error_tolerance
        hi_mass = mass + error_tolerance
        lo = self.search_binary(lo_mass)
        hi = self.search_binary(hi_mass) + 1
        return [structure for structure in self[lo:hi] if lo_mass <= structure.mass() <= hi_mass]

    def search_mass_ppm(self, mass, error_tolerance):
        """Search for the set of all items in :attr:`structures` within `error_tolerance` PPM
        of the queried `mass`.

        Parameters
        ----------
        mass : float
            The neutral mass to search for
        error_tolerance : float, optional
            The range of mass errors (in Parts-Per-Million Error) to allow

        Returns
        -------
        list
            The list of :class:`StructureRecord` instances which meet the criterion
        """
        tol = mass * error_tolerance
        return self.search_mass(mass, tol)
