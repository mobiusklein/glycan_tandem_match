from glypy import Glycan, Composition
from glypy.utils import make_struct

from ms_deisotope import CompositionListPeakDependenceGraphDeconvoluter, DistinctPatternFitter
from ms_deisotope.deconvolution import charge_range_

from brainpy import neutral_mass

from .scoring import FragmentScorer, IntensityRankScorer
from .database import Formula, CommonComposition, StructureRecord
from .preprocess import MSMSScan


class SpectrumSolution(object):
    """Describe the results of a spectrum match, including the match quality score.

    This class will convert matched peak information into a summary statistic describing
    the quality of a spectrum-structure match using information both about structural characterization
    and spectrum utilization.

    Attributes
    ----------
    assigned_peaks : list
        The list of deconvolved peaks that were assigned compositions
    matched_fragments : list
        The list of theoretical fragments for which compositions were assigned
    peaklist : ms_peak_picker.PeakIndex
        The list of experimental peaks
    score : float
        The quality of the match
    structure_record : StructureRecord
        The structure that was matched
    """
    def __init__(self, spectrum_matcher):
        self.peaklist = spectrum_matcher.peaklist
        self.structure_record = spectrum_matcher.structure_record
        self.assigned_peaks = spectrum_matcher.assigned_peaks
        self.matched_fragments = spectrum_matcher.matched_fragments
        self._structure_scorer = FragmentScorer(self.structure_record.structure, self.matched_fragments)
        self._spectrum_scorer = IntensityRankScorer(self.peaklist, self.assigned_peaks)
        self.score = self._spectrum_scorer.final_score + self._structure_scorer.final_score

    def __repr__(self):
        return "SpectrumSolution(%s, %s, %f)" % (self.peaklist, self.structure_record, self.score)


def get_fragments(assigned_peaks):
    results = set()
    for peak in assigned_peaks:
        results.update(peak.solution.base_data)
    return results


class SpectrumMatcher(object):
    """Encapsulates the process of matching experimental peaks to theoretical peaks.

    Given a :class:`ms_peak_picker.PeakIndex` from a :class:`MSMSScan` and a :class:`StructureRecord`
    instance, this type will determine the set of all compositions needed to search, handling ambiguous
    compositions and managing the searching process.

    Attributes
    ----------
    assigned_peaks : list
        The list of deconvolved peaks that were assigned compositions
    composition_list : list
        The list of compositions to attempt to assign
    matched_fragments : list
        The list of theoretical fragments for which compositions were assigned
    peaklist : ms_peak_picker.PeakIndex
        The list of experimental peaks
    structure_record : StructureRecord
        The structure to search for fragments from
    """
    def __init__(self, peaklist, structure_record):
        self.peaklist = peaklist
        self.structure_record = structure_record
        self.composition_list = None
        self.assigned_peaks = None
        self.matched_fragments = None

        self.make_composition_list()

    def make_composition_list(self):
        """Aggregate all theoretical compositions by formula so that each composition
        is present exactly once, containing all possible sources of that composition in
        a single :class:`CommonComposition` instance.

        Sets the :attr:`composition_list` attribute. Called during initialization.
        """
        self.composition_list = CommonComposition.aggregate(
            Formula(f.composition, f) for f in self.structure_record.fragment_map.values())

        # Is there anything added by including the precursor? Maybe in the graph solving step
        # f = Formula(self.structure_record.structure.total_composition(), self.structure_record)
        # self.composition_list.append(CommonComposition(f, [f.data]))

    def match(self, mass_tolerance=5e-6, charge_range=(1, 3), charge_carrier=Composition("Na").mass):
        """Perform the actual matching of peaks, fitting isotopic patterns and retrieving the matched
        structural features.

        Parameters
        ----------
        mass_tolerance : float, optional
            The parts-per-million mass error tolerance allowed. Defaults to 5ppm, entered as 5e-6
        charge_range : tuple, optional
            The a tuple defining the minimum and maximum of the range of charge states to consider
            for each composition. Defaults to positive (1, 3)
        charge_carrier : float, optional
            The mass of the charge carrying element. By default for the experimental set up described
            here, this is the mass of a sodium (Na)

        Returns
        -------
        SpectrumSolution
        """
        dec = CompositionListPeakDependenceGraphDeconvoluter(
            self.peaklist, self.composition_list,
            scorer=DistinctPatternFitter())
        self.assigned_peaks = dec.deconvolute(
            charge_range=charge_range, charge_carrier=charge_carrier,
            error_tolerance=mass_tolerance)
        self.matched_fragments = list(get_fragments(self.assigned_peaks))

        return SpectrumSolution(self)


def match(peaklist, structure_record, mass_tolerance=5e-6, charge_range=(1, 3), charge_carrier=Composition("Na").mass):
    """Matches peaks from `peaklist` to compositions derived from theoretical product ions from `structure_record`.

    Parameters
    ----------
    peaklist : ms_peak_picker.PeakIndex
        The experimental peaks to search
    structure_record : StructureRecord
        The structure to search for product ions from
    mass_tolerance : float, optional
        The mass error tolerance in ppm. Defaults to 5ppm, or 5e-6
    charge_range : tuple, optional
        The a tuple defining the minimum and maximum of the range of charge states to consider
        for each composition. Defaults to positive (1, 3)
    charge_carrier : float, optional
        The mass of the charge carrying element. By default for the experimental set up described
        here, this is the mass of a sodium (Na)

    Returns
    -------
    SpectrumSolution
    """
    return SpectrumMatcher(peaklist, structure_record).match(
        mass_tolerance=mass_tolerance, charge_range=charge_range, charge_carrier=charge_carrier)


class SolutionSet(object):
    def __init__(self, scan, solutions, name=None):
        self.scan = scan
        self.name = name
        self.solutions = solutions

    def __repr__(self):
        return "SolutionSet(%s, %r)" % (self.name, sorted(self.solutions, key=lambda x: x.score)[-1])

    def __getitem__(self, name):
        for sol in self.solutions:
            if sol.structure_record.structure.name == name:
                return sol
        raise KeyError(name)


class GlycanMSMSDatabaseSearch(object):
    def __init__(self, spectrum_dataset, structure_database, precursor_mass_tolerance=2., product_ion_mass_tolerance=5e-6,
                 maximum_charge=3, charge_carrier=Composition("Na").mass):
        self.spectrum_dataset = spectrum_dataset
        self.structure_database = structure_database
        self.precursor_mass_tolerance = precursor_mass_tolerance
        self.product_ion_mass_tolerance = product_ion_mass_tolerance
        self.maximum_charge = maximum_charge
        self.charge_carrier = charge_carrier

    def search(self):
        for spectrum in self.spectrum_dataset:
            if spectrum.has_precursor_charge():
                charge_list = [spectrum.precursor_charge]
                min_charge = 1 if spectrum.precursor_charge > 0 else -1
            else:
                if self.maximum_charge < 0:
                    charge_list = list(charge_range_(-1, self.maximum_charge))
                    min_charge = -1
                else:
                    charge_list = list(charge_range_(1, self.maximum_charge))
                    min_charge = 1

            precursor_candidates_masses = [
                (neutral_mass(spectrum.precursor_mz, c, self.charge_carrier), c) for c in charge_list
            ]

            precursor_candidate_structures = []
            for candidate_mass, candidate_charge in precursor_candidates_masses:
                precursor_candidate_structures.extend(
                    (s, candidate_charge) for s in self.structure_database.search_mass(
                        candidate_mass, self.precursor_mass_tolerance))

            solutions = []
            for precursor, charge in precursor_candidate_structures:
                solution = match(
                    spectrum.peaklist, precursor, self.product_ion_mass_tolerance,
                    charge_range=(min_charge, charge), charge_carrier=self.charge_carrier)
                solutions.append(solution)
            solution_set = SolutionSet(spectrum, solutions, spectrum.name)
            yield solution_set
