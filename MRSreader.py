"""
Read MRSolutions proprietary raw data format .MRD file into class MRSdata object
For quick check, read a particular mrd file with --input to debug

The raw .MRD file and .SPR file both have headers in the format of
':FIELDNAME ...\\r\\n' records. There are two notations for a parameter:

    pattern1: FIELDNAME value               :FOV 45
                                            :AcquisitionStartTime 13312456532890

    pattern2: FIELDNAME key, value          :OBSERVE_FREQUENCY "13C 0.0", 0.0, MHz, kHz, Hz, rx1MHz
                                            :VAR alpha, 13
                                            :SAMPLE_PERIOD sample_period, 400, 14, "25.0 KHz  40 us"

A parameter is looked up by its FIELDNAME, and by the key as well wherever the record carries one -
the sequence variables are all recorded as ':VAR', ~50 of them in one scan, and the base frequency
record is keyed by nucleus. FIELDNAME has to match in full, so 'FOV' never matches ':FOV_OFFSETS',
and of the records under one FIELDNAME the first one keyed as asked is the one read, so the key 'tr'
never picks up ':VAR tramp, 100'.
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
    DEFAULT_BASE_FREQUENCY = 74941736       # urea centered frequency in Hz
    # the binary header occupies a fixed block and the data follows it immediately
    DATA_START = 512
    # bytes one acquired point occupies, per MR Solutions data format code. Every code but 3 is
    # complex and stored as an interleaved real/imaginary pair, so the width covers both halves
    BYTES_PER_POINT = {3: 2, 16: 2, 17: 2, 18: 4, 19: 4, 20: 8, 21: 8, 22: 16}
    # dtype each format is read as. Code 3 is the only real one, the rest are read interleaved
    POINT_DTYPE = {3: 'int16', 16: 'uint8', 17: 'int8', 18: 'int16', 19: 'int16',
                   20: 'int32', 21: 'float32', 22: 'float64'}

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
        self.sample_period = 0           # in 100ns
        self.acquisition_timestamp = 0  # in 100ns since some epoch
        self.flip_angle = 0
        self.naverages = 0              # one dimension of the acquisition, read as nothing more
        self.nswitch = 1                # a divisor downstream, never 0
        self.npoints_per_switch = 0
        self.FOVoffset = [0.0, 0.0, 0.0]
        self.FOVaspect = 0.0
        self.FOV = 0.0
        self.tr = 0.0                   # in ms
        self.te = 0.0                   # in ms
        self.datatype = 0               # MR Solutions data format code
        self.rawdata = None             # np.array of shape=(nsamples, nviews, nsliceviews, nslices, nechoes, nrepetitions)
        self.parameters = ''            # sequence parameters appended as header at the end of the file
        # a missing record is reported per file, which is worth reading for one file and is noise
        # when probing a whole tree, so probing silences it
        self.quiet = False

    def _warn(self, message):
        """
        Report a record the parameter block does not carry. Silenced while probing, where the same
        handful of absences repeats once per file across a whole directory tree
        """
        if not self.quiet:
            print(message, file=sys.stderr)

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
            # a sidecar that is missing, or that carries no FREQ record, leaves the default in place
            base_frequency = self._read_frequency_from_SPR(filepath)
            if not base_frequency:
                print(f"   no base frequency read, keeping base_frequency={self.base_frequency}Hz",
                      file=sys.stderr)
                return
            self.base_frequency = base_frequency
            if DEBUG_MRSREADER:
                print(f"   setting base frequency to {self.base_frequency}Hz", file=sys.stderr)

    def parse_from_buffer(self, fdbytes):
        """
        Parse the bytes of an MR Solutions .MRD file. Split out from read_from_file(filepath) so that
        callers holding the bytes already (e.g. a member of a tar stream) can parse without a
        filesystem. When base_frequency_in_SPR comes back set, the frequency was not in these bytes and
        such a caller finishes with base_frequency = MRSdata.parse_spr(sprbytes), as read_from_file does
        Args:
            - fdbytes: complete contents of one .MRD file
        """
        try:
            shape, dend = self._read_dimensions(fdbytes)
        except ValueError as e:
            print(e, file=sys.stderr)
            return
        if self.datatype == 3:                                          # real int16, no imaginary part
            rawdata = np.frombuffer(fdbytes[self.DATA_START:dend], dtype=self.POINT_DTYPE[self.datatype])
        else:
            interleaved = np.frombuffer(fdbytes[self.DATA_START:dend], dtype=self.POINT_DTYPE[self.datatype])
            rawdata = interleaved[::2] + 1j * interleaved[1::2]
        self.rawdata = np.reshape(rawdata, shape, order='F')
        if DEBUG_MRSREADER:
            print(f"Reading {self.rawdata.shape} nsamples x nviews x nsliceviews x nslices x nechoes x nrepetitions", file=sys.stderr)
        # parameters describes settings, appended to the end of the file
        self.parameters = str(fdbytes[dend:])
        if DEBUG_MRSREADER:
            print(f"header={self.parameters}", file=sys.stderr)
        self._set_parameters()

    def _read_dimensions(self, fdbytes):
        """
        Read the fixed size binary header: the six array dimensions and the data format code. Split
        out of parse_from_buffer so that probing a file for its parameters can find where the data
        block ends without decoding the block itself
        Args:
            - fdbytes: at least the first DATA_START bytes of one .MRD file
        Returns:
            - (shape, dend): the shape rawdata takes, and the offset one past the end of the data,
              which is where the appended parameter block starts
        Raises:
            - ValueError on a data format code this reader does not know
        """
        # int(), because these are written straight into the MRD header and its serializer takes
        # Python integers only: a numpy scalar is rejected as 'not an unsigned 32-bit integer'
        self.nsamples = int(np.frombuffer(fdbytes[0:4], dtype='int32')[0])
        self.nviews = int(np.frombuffer(fdbytes[4:8], dtype='int32')[0])
        self.nsliceviews = int(np.frombuffer(fdbytes[8:12], dtype='int32')[0])
        self.nslices = int(np.frombuffer(fdbytes[12:16], dtype='int32')[0])
        self.datatype = int(np.frombuffer(fdbytes[18:20], dtype='int16')[0])
        self.nechoes = int(np.frombuffer(fdbytes[152:156], dtype='int32')[0])
        self.nrepetitions = int(np.frombuffer(fdbytes[156:160], dtype='int32')[0])
        shape = (self.nsamples, self.nviews, self.nsliceviews, self.nslices, self.nechoes, self.nrepetitions)
        if self.datatype not in self.BYTES_PER_POINT:
            raise ValueError(f"Unknown data format {self.datatype}")
        totalpts = int(np.prod(shape))
        return shape, self.DATA_START + totalpts * self.BYTES_PER_POINT[self.datatype]

    def probe_from_buffer(self, fdbytes, quiet=True):
        """
        Read the parameters of a .MRD file without decoding its data, leaving rawdata None. Grouping
        files by experiment needs the sequence name and little else, and decoding every file in a
        tree to reach a parameter block that sits after the data is what makes that expensive
        Args:
            - fdbytes: the binary header and the parameter block. The data block in between may be
              anything, including absent, as long as the parameter block starts at the same offset
            - quiet: suppress the per record 'not found' reporting, which repeats once per file
        """
        self.quiet = quiet
        try:
            _, dend = self._read_dimensions(fdbytes)
        except ValueError as e:
            self._warn(str(e))
            return
        self.parameters = str(fdbytes[dend:])
        self._set_parameters()

    def probe_from_file(self, filepath, quiet=True):
        """
        Read the parameters of a .MRD file from disk without reading its data block. Reads the fixed
        binary header, computes where the data ends from it, then seeks straight to the parameter
        block, so the cost is a few KB per file rather than the whole array.

        Unlike read_from_file this does not chase the .SPR sidecar: the base frequency plays no part
        in deciding which files belong together, and finding a sidecar costs a directory listing per
        file. A caller that wants it finishes with set_base_frequency
        Args:
            - filepath: path to the .MRD file to probe
            - quiet: suppress the per record 'not found' reporting, which repeats once per file
        """
        self.quiet = quiet
        filepath = Path(filepath)
        try:
            with open(filepath, 'rb') as fd:
                head = fd.read(self.DATA_START)
                if len(head) < self.DATA_START:
                    raise ValueError(f"file holds {len(head)} bytes, too short for a .MRD header")
                _, dend = self._read_dimensions(head)
                fd.seek(dend)
                self.parameters = str(fd.read())
        except (OSError, ValueError) as e:
            print(f"{e} at {filepath}", file=sys.stderr)
            return
        self._set_parameters()

    def set_base_frequency(self, base_frequency):
        """
        Apply a base frequency read from a .SPR sidecar elsewhere, for callers that parsed from
        bytes and so could not look the sidecar up themselves. A frequency the header already
        carried is left alone, and so is the default when the sidecar had no FREQ record
        Args:
            - base_frequency: frequency in Hz, or 0 when none was found
        """
        if not self.base_frequency_in_SPR or not base_frequency:
            return
        self.base_frequency = base_frequency
        if DEBUG_MRSREADER:
            print(f"   setting base frequency to {self.base_frequency}Hz", file=sys.stderr)

    @staticmethod
    def _parse_parameters(text, fieldname, key=None):
        """
        Acquire the values recorded for fieldname, in the two notations described at the top of the
        file. fieldname has to match in full, opened by ':' or a space and closed by a space or a ','
        Args:
            - text: repr of a parameter block, i.e. str(bytes)
            - fieldname: name of the record, e.g. 'FOV' or 'VAR'
            - key: pattern2 only, the name the record carries ahead of its value, e.g. 'alpha' or
                   '13C 0.0'. Quotes around it in the file are ignored. One fieldname can be recorded
                   many times, ':VAR' ~50 times in a scan, so every record under it is searched and
                   the first one keyed like this is the one returned
        Returns:
            - the record's comma separated values, as a list of stripped strings, the key dropped
        Raises:
            - AttributeError when fieldname is absent
            - KeyError when no record under fieldname carries key
        """
        matches = list(re.finditer(r'(?:^|[: ])' + re.escape(fieldname) + r'(?=[ ,])', text))
        if not matches:
            raise AttributeError(f"{fieldname} not found")
        for match in matches:
            # a record runs to the line ending that starts the next one, so it can span a line break
            # as ':FOV_OFFSETS 1\r\n, 0, -3.49875, 0' does, and ends at the closing quote of the repr
            tail = text[match.end():]
            end = re.search(r"\\r\\n(?=:|'|\Z)", tail)
            record = (tail[:end.start()] if end else tail).replace('\\r\\n', '')
            values = [value.strip() for value in record.split(',')]
            if key is None:
                # pattern1, the value follows the field name directly and there is nothing to choose
                # between, so the first record is the record
                return [value for value in values if value] or ['']
            if values[0].strip('"') == key:      # pattern2, the key is quoted in the file
                return values[1:]
        raise KeyError(f"no {fieldname} record is keyed {key}")

    @staticmethod
    def _parse_parameter(text, fieldname, key=None):
        """
        First of the values recorded for fieldname, which is the only one most records carry
        """
        return MRSdata._parse_parameters(text, fieldname, key)[0]

    def _set_parameters(self):
        """
        Extract value from header based on FIELDNAME notated as two patterns below:
        pattern1: FIELDNAME value               :FOV 45
                                                :AcquisitionStartTime 13312456532890

        pattern2: FIELDNAME key, value          :OBSERVE_FREQUENCY "13C 0.0", 0.0, MHz, kHz, Hz, rx1MHz
                                                :VAR alpha, 13
                                                :SAMPLE_PERIOD sample_period, 400, 14, "25.0 KHz  40 us"

        If no values were found, just sets default values on __init__
        """
        # sequence_name are stored with FIELDNAME 'SEQUENCE' for EVO2 and 'PPL' for EVO1
        # example: ':PPL C:\\smis\\dev\\Seq\\epsigre43_FB_13C.ppl' becomes 'epsigre43_FB_13C'
        for fieldname in ('SEQUENCE', 'PPL'):
            try:
                value = self._parse_parameter(self.parameters, fieldname)
            except AttributeError:
                continue                    # not this generation's field name, try the other
            else:
                # the parameter block is a repr, so a Windows path arrives with doubled backslashes
                path = value.strip('"\'').replace('\\\\', '/').replace('\\', '/')
                self.sequence_name = Path(path).stem
                if DEBUG_MRSREADER:
                    print(f"   setting sequence name to {self.sequence_name}", file=sys.stderr)
                break

        # base_frequency are stored with FIELDNAME='FREQUENCY' for EVO2 and 'OBSERVE_FREQUENCY' for EVO1
        # base_frequency = system_frequency - frequency_offset
        #
        # The key says where the frequency itself lives
            # example: ':OBSERVE_FREQUENCY "13C", 29058858.0' is an offset in Hz below the system frequency
            # example: ':FREQUENCY "13C 0.0", 0.0, frequency_base' has no offset, the frequency is in the
            #          .SPR sidecar instead, so flag it for read_from_file to pick up
        # If key="13C", the C:\smis\smis.ini file shows the system frequency via PTSMASK=190
        SYSTEM_FREQUENCY = 104000000 # smis.ini PTSmask=190 104MHz for key=13C
        for fieldname in ('FREQUENCY', 'OBSERVE_FREQUENCY'):
            try:
                offset = float(self._parse_parameter(self.parameters, fieldname, "13C"))
            except AttributeError:
                continue                    # not this generation's field name, try the other
            except KeyError:
                pass                        # the field is here but not keyed with "13C" 
            else:
                self.base_frequency = SYSTEM_FREQUENCY - int(offset)
                if DEBUG_MRSREADER:
                    print(f"   setting base frequency to {self.base_frequency}Hz", file=sys.stderr)
                break

            # If key is "13C 0.0", the offset frequency is 0
            # Get the actual frequency stored in .SPR file in the same directory
            try:
                self._parse_parameter(self.parameters, fieldname, "13C 0.0")
            except KeyError as e:
                # KeyError renders its message quoted, so unwrap it to keep the log readable
                self._warn(f"   {e.args[0]}, keeping base_frequency={self.base_frequency}Hz")
            else:
                self.base_frequency_in_SPR = True
                if DEBUG_MRSREADER:
                    print(f"   {fieldname} defers the base frequency to the .SPR sidecar", file=sys.stderr)
                break
        
        # every remaining record carries a single number, so they only differ in the attribute they
        # land on and how it is converted: (attr, fieldname, key, dtype, scale). key is None on the
        # records that hold their value directly, i.e. pattern1
        FIELDS = (
            ('naverages', 'NO_AVERAGES', 'no_averages', int, 1),                # ':NO_AVERAGES no_averages, 1'
            ('sample_period', 'SAMPLE_PERIOD', 'sample_period', int, 1),        # ':SAMPLE_PERIOD sample_period, 400, 14, "25.0 KHz 40 \xb5s"'
            ('flip_angle', 'VAR', 'alpha', int, 1),                             # ':VAR alpha, 13'
            ('tr', 'VAR', 'tr', float, 1),                                      # ':VAR tr, 60'
            ('te', 'VAR', 'te', float, 1),                                      # ':VAR te, 0'
            ('nswitch', 'VAR', 'no_switches', int, 1),                          # ':VAR no_switches, 64'
            ('npoints_per_switch', 'VAR', 'no_pts_switch', int, 1),             # ':VAR no_pts_switch, 12'
            ('FOVaspect', 'VAR', 'aspect_ratio', float, 1),                     # ':VAR aspect_ratio, 1'
            ('FOV', 'FOV', None, float, 1.0E-3),                                # ':FOV 45', in mm, stored in m
            ('acquisition_timestamp', 'AcquisitionStartTime', None, int, 1),    # ':AcquisitionStartTime 13312456532890'
        )
        for attr, fieldname, key, dtype, scale in FIELDS:
            # _parse_parameters raises on a record the block does not carry, and dtype on one it cannot
            # read, so both are handled here and the attribute keeps the value it already has
            try:
                value = dtype(self._parse_parameter(self.parameters, fieldname, key)) * scale
            except (AttributeError, KeyError, ValueError) as e:
                # KeyError renders its message quoted, so unwrap it to keep the log readable
                reason = e.args[0] if isinstance(e, KeyError) else e
                self._warn(f"   {reason}, keeping {attr}={getattr(self, attr)}")
                continue
            setattr(self, attr, value)
            if DEBUG_MRSREADER:
                print(f"   setting {attr} to {value}", file=sys.stderr)

        # FOV offsets are the one record holding more than one value, ':FOV_OFFSETS 1\r\n, 0, -3.49875,
        # 0' - a count followed by the three offsets, in mm in the file and stored in m
        try:
            values = self._parse_parameters(self.parameters, 'FOV_OFFSETS')
            offsets = [float(value) / 1000.0 for value in values[1:4]]
            if len(offsets) != 3:
                raise ValueError(f"FOV_OFFSETS holds {len(offsets)} offsets, not 3")
        except (AttributeError, ValueError) as e:
            self._warn(f"   {e}, keeping FOVoffset={self.FOVoffset}")
        else:
            self.FOVoffset = offsets
            if DEBUG_MRSREADER:
                print(f"   setting FOV offsets (m) to {self.FOVoffset}", file=sys.stderr)

    # ---------- finding the EPSI readout sampling window ----------------------

    def switch_layout(self):
        """
        How one EPSI readout divides into gradient switches.

        no_switches and no_pts_switch are recorded separately from the sample count, so the
        acquisition holds more points per switch than the sequence calls usable: the extra ones
        are the gradient ramps either side of the flat top.
        Returns:
            - (nswitch, total points in one switch, points the sequence keeps per switch)
        """
        nswitch = max(self.nswitch, 1)
        total = self.nsamples // nswitch
        return nswitch, total, (self.npoints_per_switch or total)

    def switch_profile(self, spectral=True, view=None, repeat=None):
        """
        How much signal sits at each position within a readout switch.

        Every sample of the readout falls at one of `total` positions inside a gradient switch,
        and the gradient echo, the point where kx crosses zero, sits at one of them. Summing raw
        magnitude over the switches locates it only when the object fills the field of view. For a
        compact object |k-space| is close to flat and the echo hides in the phase instead, so the
        default first transforms along the switch axis and keeps the strongest spectral line,
        which is a coherent sum over all the switches and lifts the echo well clear of the noise.
        Args:
            - spectral: aggregate through the spectral transform rather than by raw magnitude
            - view: phase encode line to read, or None for the one carrying the most signal
            - repeat: index into the folded slice/echo/repetition axis, or None to sum over it
        Returns:
            - (profile of length total, peak over median of the profile)
        Raises:
            - ValueError when called on an object whose data block was never read
        """
        if self.rawdata is None:
            raise ValueError("no data to profile; read_from_file rather than probe_from_file")
        nswitch, total, _ = self.switch_layout()
        # (switch, position in switch, view, everything else). The sample axis comes first, so
        # splitting it in C order is exactly the switch-major order the samples were acquired in.
        # The slice, sliceview and echo axes are single valued in an EPSI scan and are folded in
        # with the repetitions, since for this purpose they are all just repeats
        cube = self.rawdata[:nswitch * total].reshape(nswitch, total, self.nviews, -1)

        if spectral:
            spectrum = np.fft.fft(cube, axis=0)
            power = (np.abs(spectrum) ** 2).sum(axis=(1, 2, 3))
            band = np.abs(spectrum[int(np.argmax(power))])          # (position, view, repeat)
        else:
            band = np.abs(cube).sum(axis=0)                         # (position, view, repeat)

        if view is None:
            view = int(np.argmax(band.sum(axis=(0, 2))))
        lines = band[:, view, :]
        profile = lines[:, repeat] if repeat is not None else lines.sum(axis=1)
        median = np.median(profile)
        return profile, (float(profile.max() / median) if median else np.inf)

    def sliding_window(self, profile, kept):
        """
        Score every candidate sampling window of `kept` consecutive positions.

        Windows wrap, because the last position of one switch is followed by the first position of
        the next, so a window may legitimately straddle the boundary.
        Args:
            - profile: signal per position within the switch, from switch_profile
            - kept: width of the window, i.e. the points the reconstruction keeps per switch
        Returns:
            - (signal summed inside each window, position of the profile peak within each window),
              both indexed by the window's first position
        """
        total = len(profile)
        offsets = np.arange(kept)
        signal = np.array([profile[(start + offsets) % total].sum() for start in range(total)])
        # where the echo lands inside the window; the window is right when this is its middle
        peak_at = (int(np.argmax(profile)) - np.arange(total)) % total
        return signal, peak_at

    def window_report(self, spectral=True, view=None, repeat=None):
        """
        Pick the sampling window, both ways, and translate it for the reconstruction.

        Two criteria, because they can disagree and the disagreement is the interesting part: the
        window holding the most signal, and the window whose middle sits on the echo. mrd2recon
        addresses the same window as an offset from discard_pre, so that conversion is reported
        alongside, as the EPSIGRE_LEADING_PAD that would select it.
        Returns:
            - dict of the profile, the two winning window starts and the pads they imply
        """
        nswitch, total, kept = self.switch_layout()
        profile, snr = self.switch_profile(spectral=spectral, view=view, repeat=repeat)
        signal, peak_at = self.sliding_window(profile, kept)
        middle = (kept - 1) / 2

        by_signal = int(np.argmax(signal))
        # distance from the middle measured the short way round, so a wrapped window is not
        # penalised for having its peak reported as position 27 rather than -1
        from_middle = np.abs((peak_at - middle + total / 2) % total - total / 2)
        by_centre = int(np.argmin(from_middle))
        discard_pre = (total - kept) // 2

        def pad_for(start):
            """
            The EPSIGRE_LEADING_PAD that makes mrd2recon read this window.

            It addresses a window as discard_pre - pad, so the pad is that difference, wrapped
            the short way round the switch: the positions are cyclic, so a window starting at 29
            of 34 is 5 before the boundary rather than 29 after it
            """
            return int((discard_pre - start + total // 2) % total - total // 2)

        return dict(profile=profile, signal=signal, peak_at=peak_at, snr=snr,
                    nswitch=nswitch, total=total, kept=kept, discard_pre=discard_pre,
                    peak=int(np.argmax(profile)),
                    by_signal=by_signal, by_centre=by_centre,
                    pad_by_signal=pad_for(by_signal),
                    pad_by_centre=pad_for(by_centre))

    def plot_switch_profile(self, spectral=True, view=None, repeat=None, savepath=None, show=True):
        """
        Plot the sliding window scan over the positions within a readout switch.

        Three panels: the signal at each position with the chosen windows drawn on it, the sliding
        window scan itself, and the same profile per switch so that an echo whose position drifts
        along the echo train shows up rather than being averaged away.
        Args:
            - savepath: write the figure here instead of, or as well as, showing it
            - show: open a window, which blocks until it is closed
        Returns:
            - the window_report dict, so a caller can act on the numbers it plotted
        """
        import matplotlib.pyplot as plt

        report = self.window_report(spectral=spectral, view=view, repeat=repeat)
        total, kept = report['total'], report['kept']
        profile, signal = report['profile'], report['signal']
        middle = (kept - 1) / 2
        positions = np.arange(total)

        figure, axes = plt.subplots(3, 1, figsize=(11, 10))
        title = (f"{self.sequence_name or 'unknown sequence'}: {report['nswitch']} switches of "
                 f"{total} points, {kept} kept, peak/median {report['snr']:.2f}")
        figure.suptitle(title)

        axes[0].plot(positions, profile, 'o-', color='C0')
        axes[0].axvline(report['peak'], color='C1', lw=2, label=f"echo peak at {report['peak']}")
        for start, colour, name in ((report['by_signal'], 'C2', 'most signal'),
                                    (report['by_centre'], 'C3', 'echo centred')):
            # drawn as the positions it covers, so a window that wraps appears at both ends
            covered = (start + np.arange(kept)) % total
            axes[0].plot(covered, profile[covered], 'o', ms=11, mfc='none', color=colour,
                         label=f"{name}: start {start}, pad "
                               f"{report['discard_pre'] - start}")
        axes[0].set_xlabel(f"position within the {total} point switch")
        axes[0].set_ylabel("signal")
        axes[0].legend(fontsize=8)

        axes[1].plot(positions, signal / signal.max(), 'o-', color='C0', label='signal in window')
        axes[1].plot(positions, report['peak_at'] / total, 's--', ms=3, color='C4',
                     label='where the peak lands in the window')
        axes[1].axhline(middle / total, color='C3', ls=':', label='window middle')
        axes[1].axvline(report['by_signal'], color='C2', lw=2)
        axes[1].axvline(report['by_centre'], color='C3', lw=2)
        axes[1].set_xlabel("first position of the window")
        axes[1].set_ylabel("normalised")
        axes[1].legend(fontsize=8)

        nswitch = report['nswitch']
        cube = self.rawdata[:nswitch * total].reshape(nswitch, total, self.nviews, -1)
        axes[2].imshow(np.abs(cube).sum(axis=(2, 3)), aspect='auto', origin='lower',
                       interpolation='nearest')
        axes[2].axvline(report['by_centre'], color='C3', lw=1.5)
        axes[2].axvline((report['by_centre'] + kept - 1) % total, color='C3', lw=1.5)
        axes[2].set_xlabel(f"position within the switch")
        axes[2].set_ylabel("switch")

        figure.tight_layout()
        if savepath:
            figure.savefig(savepath, dpi=110)
            print(f"wrote {savepath}", file=sys.stderr)
        if show:
            plt.show()
        return report

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

if __name__ == '__main__':
    DEBUG_MRSREADER = True
    parser = argparse.ArgumentParser(description='Read an MRS .MRD file and print the parsed parameters')
    parser.add_argument('-i', '--input', type=Path, required=True,
                        help='path to input .MRD file containing MRS data')
    parser.add_argument('-w', '--window', action='store_true',
                        help='scan the sampling window within an EPSI readout switch and plot it')
    parser.add_argument('--magnitude', action='store_true',
                        help='profile by raw magnitude instead of the spectral transform, which '
                             'only finds the echo when the object fills the field of view')
    parser.add_argument('--view', type=int, default=None,
                        help='phase encode line to profile (default: the brightest)')
    parser.add_argument('--repeat', type=int, default=None,
                        help='single repetition to profile (default: sum over all of them)')
    parser.add_argument('--save', type=Path, default=None, help='write the figure to this path')
    parser.add_argument('--no-show', action='store_true', help='do not open a plot window')
    args = parser.parse_args()
    mrs = MRSdata()
    mrs.read_from_file(args.input)
    if args.window:
        report = mrs.plot_switch_profile(spectral=not args.magnitude, view=args.view,
                                         repeat=args.repeat, savepath=args.save,
                                         show=not args.no_show)
        print(f"{report['nswitch']} switches of {report['total']} points, keeping {report['kept']}, "
              f"discard_pre={report['discard_pre']}, profile peak/median {report['snr']:.2f}")
        print(f"echo peak at position {report['peak']}")
        print(f"most signal:  window starts at {report['by_signal']}, echo lands at "
              f"{report['peak_at'][report['by_signal']]} of {report['kept']}, "
              f"EPSIGRE_LEADING_PAD = {report['pad_by_signal']}")
        print(f"echo centred: window starts at {report['by_centre']}, echo lands at "
              f"{report['peak_at'][report['by_centre']]} of {report['kept']}, "
              f"EPSIGRE_LEADING_PAD = {report['pad_by_centre']}")
        if report['snr'] < 1.5:
            print("the profile is nearly flat, so this scan carries too little signal to place "
                  "the window; check it against one that does", file=sys.stderr)
