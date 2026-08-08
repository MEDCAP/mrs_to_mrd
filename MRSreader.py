"""
Read MRSolutions proprietary raw data format .MRD file into class MRSdata object
For quick check, read a particular mrd file with --input to debug

The settings appended after the raw data, and the .SPR sidecar next to it, are both a sequence of
':FIELDNAME ...\\r\\n' records. There are two notations for a parameter:

    pattern1: FIELDNAME value               :FOV 45
                                            :AcquisitionStartTime 13312456532890

    pattern2: FIELDNAME key, value          :OBSERVE_FREQUENCY "13C 0.0", 0.0, MHz, kHz, Hz, rx1MHz
                                            :VAR alpha, 13
                                            :SAMPLE_PERIOD sample_period, 400, 14, "25.0 KHz  40 us"

A parameter is looked up by the name that sits immediately before its value, so pattern1 is looked up
on the field name and pattern2 on the key: 'alpha' rather than 'VAR', 'sample_period' rather than
'SAMPLE_PERIOD'. One scan holds ~50 ':VAR' records, so the field name cannot single one out. Only the
frequency needs both, because there its key says where the frequency lives rather than naming a
variable, so _parse_parameters takes a key to check against.

Either way the name has to match in full: 'FOV' matches ':FOV 45' but never ':FOV_OFFSETS 1' or
':VAR FOVf, 12', and 'tr' matches ':VAR tr, 60' but never ':VAR tramp, 100'.
"""

from pathlib import Path
import numpy as np
import re
import sys
import argparse

# enable below to print the obtained variables
DEBUG_MRSREADER = False

class MRSdata:
    # used when neither the header nor a .SPR sidecar gives a frequency
    DEFAULT_BASE_FREQUENCY = 74941736       # urea centered frequency

    def __init__(self):
        self.nsamples = 0
        self.nviews = 0
        self.nsliceviews = 0
        self.nslices = 0
        self.nechoes = 0
        self.nrepetitions = 0
        self.base_frequency = self.DEFAULT_BASE_FREQUENCY
        self.base_frequency_in_SPR = False  # set when the header defers the frequency to the sidecar
        self.sequence_name = ''
        self.sampleperiod = 0           # in 100ns
        self.nslc = 0
        self.acquisition_timestamp = 0  # in 100ns since some epoch
        self.flip_angle = 0
        self.naverages = 1              # a phantom marker downstream, never 0
        self.nswitch = 1                # a divisor downstream, never 0
        self.npoints_per_switch = 0
        self.FOVoffset = [0.0, 0.0, 0.0]
        self.FOVaspect = 0.0
        self.FOV = 0.0
        self.tr = 0.0                   # in ms
        self.datatype = 0               # MR Solutions data format code
        self.rawdata = None             # (nsamples, nviews, nsliceviews, nslices, nechoes, nrepetitions)
        self.parameters = ''            # sequence settings appended to the end of the file

    def read_from_file(self, filepath):
        """
        Read the MRS file into class fields by passing the filepath. The base frequency comes from a
        .SPR sidecar in the same folder when the header defers to it
        """
        filepath = Path(filepath)
        try:
            if DEBUG_MRSREADER:
                print(f"Reading MRS file at {filepath}", file=sys.stderr)
            self.parse_from_buffer(filepath.read_bytes())
        except Exception as e:
            print(f"{e} at {filepath}", file=sys.stderr)
            return
        if self.base_frequency_in_SPR:
            self.set_base_frequency(self._read_frequency_from_SPR(filepath))

    def parse_from_buffer(self, fdbytes):
        """
        Parse the bytes of an MR Solutions .MRD file. Split out from read_from_file(filepath) so that
        callers holding the bytes already (e.g. a member of a tar stream) can parse without a
        filesystem. When base_frequency_in_SPR comes back set, the frequency was not in these bytes
        and such a caller finishes with set_base_frequency(MRSdata.parse_spr(sprbytes))
        Args:
            - fdbytes: complete contents of one .MRD file
        """
        # kept as python ints rather than the numpy scalars frombuffer returns, since these are
        # loop bounds and are written straight into MRD header fields, whose serializer takes an
        # int or a same-width numpy scalar but rejects the int32 read here
        self.nsamples = int(np.frombuffer(fdbytes[0:4], dtype='int32')[0])
        self.nviews = int(np.frombuffer(fdbytes[4:8], dtype='int32')[0])
        self.nsliceviews = int(np.frombuffer(fdbytes[8:12], dtype='int32')[0])
        self.nslices = int(np.frombuffer(fdbytes[12:16], dtype='int32')[0])
        self.datatype = int(np.frombuffer(fdbytes[18:20], dtype='int16')[0])
        self.nechoes = int(np.frombuffer(fdbytes[152:156], dtype='int32')[0])
        self.nrepetitions = int(np.frombuffer(fdbytes[156:160], dtype='int32')[0])
        shape = (self.nsamples, self.nviews, self.nsliceviews, self.nslices, self.nechoes, self.nrepetitions)
        totalpts = int(np.prod(shape))
        dstart = 512
        if self.datatype == 3:                                          # real int16, no imaginary part
            dend = dstart + totalpts * 2
            rawdata = np.frombuffer(fdbytes[dstart:dend], dtype='int16')
        elif self.datatype == 16:                                       # complex uint8
            dend = dstart + totalpts * 2
            interleaved = np.frombuffer(fdbytes[dstart:dend], dtype='uint8')
            rawdata = interleaved[::2] + 1j * interleaved[1::2]
        elif self.datatype == 17:                                       # complex int8
            dend = dstart + totalpts * 2
            interleaved = np.frombuffer(fdbytes[dstart:dend], dtype='int8')
            rawdata = interleaved[::2] + 1j * interleaved[1::2]
        elif self.datatype == 18 or self.datatype == 19:                # complex int16
            dend = dstart + totalpts * 4
            interleaved = np.frombuffer(fdbytes[dstart:dend], dtype='int16')
            rawdata = interleaved[::2] + 1j * interleaved[1::2]
        elif self.datatype == 20:                                       # complex int32
            dend = dstart + totalpts * 8
            interleaved = np.frombuffer(fdbytes[dstart:dend], dtype='int32')
            rawdata = interleaved[::2] + 1j * interleaved[1::2]
        elif self.datatype == 21:                                       # complex float32
            dend = dstart + totalpts * 8
            interleaved = np.frombuffer(fdbytes[dstart:dend], dtype='float32')
            rawdata = interleaved[::2] + 1j * interleaved[1::2]
        elif self.datatype == 22:                                       # complex float64
            dend = dstart + totalpts * 16
            interleaved = np.frombuffer(fdbytes[dstart:dend], dtype='float64')
            rawdata = interleaved[::2] + 1j * interleaved[1::2]
        else:
            print("Unknown data format", file=sys.stderr)
            return
        self.rawdata = np.reshape(rawdata, shape, order='F')
        if DEBUG_MRSREADER:
            print(f"Reading {self.rawdata.shape} nsamples x nviews x nsliceviews x nslices x nechoes x nrepetitions",
                  file=sys.stderr)
        # parameters describes settings, appended to the end of the file
        self.parameters = str(fdbytes[dend:])
        self._set_parameters()

    @staticmethod
    def _parse_parameters(text, fieldname, key=None):
        """
        Acquire the values recorded for fieldname, in the two notations described at the top of the
        file. fieldname has to match in full, opened by ':' or a space and closed by a space or a ',',
        and every name looked up here appears once, so there is nothing to choose between
        Args:
            - text: repr of a parameter block, i.e. str(bytes)
            - fieldname: name immediately preceding the value, e.g. 'FOV' or 'sample_period'
            - key: pattern2 only, the key the record has to carry, e.g. '13C 0.0'. Quotes around the
                   key in the file are ignored
        Returns:
            - the record's comma separated values, as a list of stripped strings
        Raises:
            - AttributeError when fieldname is absent
            - KeyError when key is given and the record carries a different one
        """
        match = re.search(r'(?:^|[: ])' + re.escape(fieldname) + r'(?=[ ,])', text)
        if match is None:
            raise AttributeError(f"{fieldname} not found")
        # a record runs to the line ending that starts the next one, so it can span a line break as
        # ':FOV_OFFSETS 1\r\n, 0, -3.49875, 0' does, and ends at the closing quote of the repr
        tail = text[match.end():]
        end = re.search(r"\\r\\n(?=:|'|\Z)", tail)
        record = (tail[:end.start()] if end else tail).replace('\\r\\n', '')
        values = [value.strip() for value in record.split(',')]
        if key is None:
            # pattern1, the value follows the name directly. Drop the empty leading item left by the
            # ',' of ':VAR alpha, 13', where the name is what precedes the comma
            return [value for value in values if value] or ['']
        if values[0].strip('"') != key:          # pattern2, the key is quoted in the file
            raise KeyError(f"{fieldname} is keyed {values[0]}, not {key}")
        return values[1:]

    @staticmethod
    def _parse_parameter(text, fieldname, key=None):
        """
        First of the values recorded for fieldname, which is the only one most records carry
        """
        return MRSdata._parse_parameters(text, fieldname, key)[0]

    def _set_parameter(self, attr, fieldname, dtype=int, scale=1):
        """
        Look up fieldname in the parameter block and store its value in attr, converted with dtype and
        multiplied by scale. A field that is absent, or that carries something dtype cannot read,
        leaves the default from __init__ in place
        """
        try:
            value = dtype(self._parse_parameter(self.parameters, fieldname)) * scale
        except (AttributeError, ValueError) as e:
            print(f"   {e}, keeping {attr}={getattr(self, attr)}", file=sys.stderr)
            return
        setattr(self, attr, value)
        if DEBUG_MRSREADER:
            print(f"   setting {attr} to {value}", file=sys.stderr)

    def _set_parameters(self):
        """
        Extract the sequence settings from the parameter block appended after the raw data.
        Every field is looked up on its own so it can be traced back to the record it comes from, and
        anything the block does not carry keeps the default from __init__
        """
        self._set_sequence_name()
        self._set_base_frequency()
        # ':SAMPLE_PERIOD sample_period, 400, 14, "25.0 KHz  40 \xb5s"'
        self._set_parameter('sampleperiod', 'sample_period', int)
        # ':NO_SLICES no_slices, 1'
        self._set_parameter('nslc', 'no_slices', int)
        # ':NO_AVERAGES no_averages, 1'
        self._set_parameter('naverages', 'no_averages', int)
        # ':VAR alpha, 13'
        self._set_parameter('flip_angle', 'alpha', int)
        # ':VAR tr, 60'
        self._set_parameter('tr', 'tr', float)
        # ':VAR no_switches, 64'
        self._set_parameter('nswitch', 'no_switches', int)
        # ':VAR no_pts_switch, 12'
        self._set_parameter('npoints_per_switch', 'no_pts_switch', int)
        # ':VAR aspect_ratio, 1'
        self._set_parameter('FOVaspect', 'aspect_ratio', float)
        # ':FOV 45', in mm, stored in m
        self._set_parameter('FOV', 'FOV', float, scale=1.0E-3)
        # ':AcquisitionStartTime 13312456532890'
        self._set_parameter('acquisition_timestamp', 'AcquisitionStartTime', int)
        self._set_fov_offsets()

    def _set_sequence_name(self):
        """
        Set the sequence name from the sequence file path, recorded at the 'PPL' field on EVO1 and the
        'SEQUENCE' field on EVO2. Only the bare name is kept, so
        ':PPL C:\\smis\\dev\\Seq\\epsigre43_FB_13C.ppl' becomes 'epsigre43_FB_13C'
        """
        for fieldname in ('SEQUENCE', 'PPL'):
            try:
                value = self._parse_parameter(self.parameters, fieldname)
            except AttributeError:
                continue
            # the parameter block is a repr, so a Windows path arrives with doubled backslashes
            path = value.strip('"\'').replace('\\\\', '/').replace('\\', '/')
            self.sequence_name = Path(path).stem
            if DEBUG_MRSREADER:
                print(f"   setting sequence name to {self.sequence_name}", file=sys.stderr)
            return
        print(f"   SEQUENCE/PPL not found, keeping sequence_name='{self.sequence_name}'", file=sys.stderr)

    def _set_base_frequency(self):
        """
        Set the base frequency from the frequency record, named 'OBSERVE_FREQUENCY' on EVO1 and
        'FREQUENCY' on EVO2. Its key says where the frequency itself lives:
            - keyed '13C', e.g. ':FREQUENCY "13C", 1234' - an offset in Hz below the transmitter
              board's own frequency, which is what PTSMASK_190_FREQUENCY records
            - keyed '13C 0.0', e.g. ':OBSERVE_FREQUENCY "13C 0.0", 0.0, MHz' - the frequency is in the
              .SPR sidecar instead, so flag it for read_from_file to pick up
        A record keyed for another nucleus, or no record at all, leaves the default from __init__
        """
        PTSMASK_190_FREQUENCY = 104000000 # smis.ini PTSmask=190 104MHz for key=13C
        for fieldname in ('FREQUENCY', 'OBSERVE_FREQUENCY'):    # two version of fieldnames for EVO1 and EVO2
            try:
                offset = float(self._parse_parameter(self.parameters, fieldname, "13C"))
            except AttributeError:
                continue                    # not this generation's field name, try the other
            except (KeyError, ValueError):
                pass                        # the field is here but not keyed with a bare nucleus
            else:
                self.base_frequency = PTSMASK_190_FREQUENCY - int(offset)
                if DEBUG_MRSREADER:
                    print(f"   setting base frequency to {self.base_frequency}Hz", file=sys.stderr)
                return
            try:
                self._parse_parameter(self.parameters, fieldname, "13C 0.0")
            except KeyError as e:
                # KeyError renders its message quoted, so unwrap it to keep the log readable
                print(f"   {e.args[0]}, keeping base_frequency={self.base_frequency}Hz", file=sys.stderr)
                return
            self.base_frequency_in_SPR = True
            if DEBUG_MRSREADER:
                print(f"   {fieldname} defers the base frequency to the .SPR sidecar", file=sys.stderr)
            return
        print(f"   FREQUENCY/OBSERVE_FREQUENCY not found, keeping base_frequency={self.base_frequency}Hz",
              file=sys.stderr)

    def _set_fov_offsets(self):
        """
        Set the three FOV offsets (mm in the file, stored in m). Looked up on its own because the
        record holds a count before the three of them, ':FOV_OFFSETS 1\\r\\n, 0, -3.49875, 0'
        """
        try:
            values = self._parse_parameters(self.parameters, 'FOV_OFFSETS')
            offsets = [float(value) / 1000.0 for value in values[1:4]]
            if len(offsets) != 3:
                raise ValueError(f"FOV_OFFSETS holds {len(offsets)} offsets, not 3")
        except (AttributeError, ValueError) as e:
            print(f"   {e}, keeping FOVoffset={self.FOVoffset}", file=sys.stderr)
            return
        self.FOVoffset = offsets
        if DEBUG_MRSREADER:
            print(f"   setting FOV offsets (m) to {self.FOVoffset}", file=sys.stderr)

    @staticmethod
    def parse_spr(sprbytes):
        """
        Read the base frequency out of a .SPR sidecar, recorded as ':EDITTEXT FREQ 74.942486000000' in
        MHz. Takes the bytes rather than a path so that a caller holding them already (e.g. a member of
        a tar stream) can read the frequency without a filesystem
        Args:
            - sprbytes: complete contents of one .SPR file
        Returns:
            - base frequency in Hz, or 0 when the file carries no FREQ record
        """
        try:
            return int(float(MRSdata._parse_parameter(str(sprbytes), 'FREQ')) * 1.0E+6 + 0.5)
        except (AttributeError, ValueError) as e:
            print(f"   no base frequency in the .SPR sidecar: {e}", file=sys.stderr)
            return 0

    def _read_frequency_from_SPR(self, filepath):
        """
        Find file with extension .SPR in the same directory as .MRD and read its FREQ record
        Args:
            - filepath: path to the .MRD file whose sidecar to look for
        Returns:
            - base frequency in Hz, or 0 when no sidecar carries one
        """
        base_frequency = 0
        try:
            for auxfile in Path(filepath).parent.iterdir():
                if auxfile.is_file() and auxfile.suffix == '.SPR':
                    freq = self.parse_spr(auxfile.read_bytes())
                    if freq:                # keep an earlier hit if this SPR has no FREQ record
                        base_frequency = freq
        except OSError as e:
            print(f"{e} looking for a .SPR sidecar next to {filepath}", file=sys.stderr)
        return base_frequency

    def set_base_frequency(self, base_frequency=0):
        """
        Store the base frequency in Hz, keeping the default from __init__ when 0 is given
        """
        if not base_frequency:
            print(f"   no base frequency read, keeping base_frequency={self.base_frequency}Hz",
                  file=sys.stderr)
            return
        self.base_frequency = base_frequency
        if DEBUG_MRSREADER:
            print(f"   setting base frequency to {self.base_frequency}Hz", file=sys.stderr)

if __name__ == '__main__':
    DEBUG_MRSREADER = True
    parser = argparse.ArgumentParser(description='Read an MRS .MRD file and print the parsed parameters')
    parser.add_argument('-i', '--input', type=Path, required=True,
                        help='path to input .MRD file containing MRS data')
    args = parser.parse_args()
    mrs = MRSdata()
    mrs.read_from_file(args.input)
