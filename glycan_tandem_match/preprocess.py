import os

from collections import namedtuple
from pyteomics import mzml
from ms_peak_picker import pick_peaks as _pick_peaks, scan_filter
from ms_peak_picker.utils import Base
from glob import glob


class MSMSScan(Base):
    def __init__(self, name, precursor_mz, precursor_charge, peaklist):
        self.name = name
        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge
        self.peaklist = peaklist

    def has_precursor_charge(self):
        return self.precursor_charge is not None


def pick_peaks(mzml_file):
    scan = mzml.read(mzml_file).next()
    return process_mzml_scan(scan, savgol_window_length=7)


def process_mzml_scan(scan, savgol_window_length=None, remove_baseline=True):
    filters = []
    if remove_baseline:
        filters.append(scan_filter.fticr_remove_baseline)
    if savgol_window_length is not None:
        savgol_filter = scan_filter.SavitskyGolayFilter(window_length=savgol_window_length)
        filters.append(savgol_filter)
    mz_array = scan["m/z array"]
    intensity_array = scan["intensity array"]
    precursor = scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]
    precursor_mz = precursor['selected ion m/z']
    try:
        precursor_charge = precursor["charge state"]
    except:
        precursor_charge = None
    return MSMSScan(
        scan['id'],
        precursor_mz, precursor_charge,
        _pick_peaks(mz_array, intensity_array, fit_type='lorenztian', transforms=filters))


class SpectrumReader(object):
    def __init__(self, file_path):
        self.file_path = file_path
        self._iter = None
        self._make_iter()
        self.name = os.path.basename(file_path)

    def _make_iter(self):
        raise NotImplementedError()

    def next(self):
        return self._iter.next()

    def __iter__(self):
        self._make_iter()
        for x in self._iter:
            yield x

    def __repr__(self):
        return "{self.__class__.__name__}({self.name})".format(self=self)


class MzMLReader(SpectrumReader):
    def __init__(self, file_path, savgol_window_length=7, remove_baseline=True):
        self.savgol_window_length = savgol_window_length
        self.remove_baseline = remove_baseline
        super(MzMLReader, self).__init__(file_path)

    def _make_iter(self):
        self._reader = mzml.read(self.file_path)
        self._iter = (
            process_mzml_scan(
                scan, savgol_window_length=self.savgol_window_length,
                remove_baseline=self.remove_baseline)
            for scan in self._reader)
