from mepi import lib_log
import numpy as np
from casatools import msmetadata
from casatools import table as tbtool
from contextlib import contextmanager

log = lib_log.log

# this is missing in casatools.table, so we wrap it in a context manager to ensure proper closing
@contextmanager
def open_table(path, readonly=False):
    tb = tbtool()
    tb.open(path, nomodify=readonly)
    try:
        yield tb
    finally:
        tb.close()

class MS:
    def __init__(self, msfile):
        self.msfile = msfile
        self.log = log
        log.info(f"Opening MS: {msfile}")
        with open_table(self.msfile + '::SPECTRAL_WINDOW') as t:
            self.freq_center = t.getcol('CHAN_FREQ').T[0].mean() # Hz - .T for casa tables, not needed for pyrap tables
        self.band = self.get_band()

    # Known band tokens in the SPECTRAL_WINDOW NAME (case-insensitive)
    _BAND_TOKENS = ('UHF', 'S2', 'S1', 'S0', 'S', 'L')  # order matters: longer/specific first

    def get_band(self):
        """
        Identify the MeerKAT receiver band from the SPECTRAL_WINDOW NAME metadata.

        MeerKAT MSes store the receiver name in the NAME column of the
        SPECTRAL_WINDOW subtable (e.g. 'MeerKAT:UHF', 'MeerKAT:L', 'MeerKAT:S').
        Falls back to central frequency if the name is absent or unrecognised.

        Returns
        -------
        str
            Band name: one of 'UHF', 'L', 'S0', 'S1', 'S2'.
            Returns 'UNKNOWN' if detection fails.
        """
        with open_table(self.msfile + '::SPECTRAL_WINDOW') as t:
            spw_name = t.getcol('NAME')[0].upper()

        for token in self._BAND_TOKENS:
            if token in spw_name:
                self.log.info(f"Band detection: SPW NAME '{spw_name}' -> {token}")
                return token

        # Fallback: use central frequency
        self.log.warning(f"Could not parse band from SPW NAME '{spw_name}', falling back to frequency.")
        if self.freq_center < 1.0e9:
            band = 'UHF'
        elif self.freq_center < 1.8e9:
            band = 'L'
        elif self.freq_center < 2.6e9:
            band = 'S0'
        elif self.freq_center < 3.1e9:
            band = 'S1'
        else:
            band = 'S2'
        self.log.info(f"Band detection: central freq {self.freq_center*1e-6:.0f} MHz -> {band}")
        return band

    def get_good_spw(self):
        """
        Return the SPW selection string for RFI-clean channels of the current band.

        Returns
        -------
        str
            CASA-style SPW selection (e.g. '0:842~869MHz').
        """
        good_chans = {
            'UHF': '0:842~869MHz',
            'L':   '0:1326~1367MHz',
            'S0':  '0:2010~2096MHz',
            'S1':  '0:2517~2603MHz',
            'S2':  '0:2876~2962MHz',
        }
        spw_good = good_chans.get(self.band)
        if spw_good is None:
            raise ValueError(f"No good channel selection defined for band '{self.band}'")
        self.log.info(f"Good channels for band {self.band}: {spw_good}")
        return spw_good

    def _find_fields_for_intent(self, intent_pattern):
        """
        Return a list of field names whose scan intent matches intent_pattern.
        Returns a single name string if only one match, otherwise a list.
        """
        msmd = msmetadata()
        msmd.open(self.msfile)
        try:
            fields = list(msmd.fieldsforintent(intent_pattern, asnames=True))
        finally:
            msmd.close()
        if not fields:
            self.log.warning(f"No field found for intent '{intent_pattern}'")
            return None
        self.log.info(f"Found fields for intent '{intent_pattern}': {fields}")
        return fields

    def find_bandpass_calibrator(self):
        """Return the field name used as bandpass calibrator (intent CALIBRATE_BANDPASS*)."""
        return self._find_fields_for_intent("CALIBRATE_BANDPASS*")

    def find_phase_calibrator(self):
        """Return the field name used as phase calibrator (intent CALIBRATE_PHASE*)."""
        return self._find_fields_for_intent("CALIBRATE_PHASE*")

    def find_targets(self):
        """Return a list of field names used as target (intent TARGET*)."""
        return self._find_fields_for_intent("TARGET*")

    # Known 3C286 / J1331+3030 aliases used as pol calibrator at MeerKAT
    _POL_CAL_NAMES = ('J1331+3030', '3C286', '3c286')

    def find_pol_calibrator(self):
        """Return the field name used as polarisation calibrator.

        First checks if any field named J1331+3030 (or its aliases) is present
        in the MS, since MeerKAT observations sometimes lack a CALIBRATE_POL
        intent even when 3C286 is observed. Falls back to intent-based lookup.
        """
        msmd = msmetadata()
        msmd.open(self.msfile)
        try:
            all_fields = list(msmd.fieldnames())
        finally:
            msmd.close()

        for name in self._POL_CAL_NAMES:
            if name in all_fields:
                self.log.info(f"Polarisation calibrator found by name: {name}")
                return [name]

        return self._find_fields_for_intent("CALIBRATE_POL*")

    def get_spw_noedges(self):
        with open_table(self.msfile+"::SPECTRAL_WINDOW") as t:
            chan_freqs = t.getcol("CHAN_FREQ")
            nchan = chan_freqs.size
            spw_selection = f"0:{nchan//20}~{nchan - nchan//20 - 1}" # select all channels but the first and last 5%
            self.log.debug(f"SPW selection (no edges): {spw_selection} - total {nchan} channels, excluding first and last {nchan//20} channels")
            return spw_selection

    def find_reference_antenna(self):
        """
        Return the name of the most distant unflagged antenna from the array centre.

        Reads POSITION (ITRF XYZ) and FLAG_ROW from the ANTENNA subtable,
        computes each unflagged antenna's distance from the centroid of all
        unflagged positions, and returns the one furthest away.
        """
        with open_table(self.msfile + '::ANTENNA') as t:
            names    = t.getcol('NAME')
            positions = t.getcol('POSITION')   # shape (nant, 3), ITRF XYZ in metres
            flags    = t.getcol('FLAG_ROW')

        unflagged = [(name, pos) for name, pos, flag in zip(names, positions, flags) if not flag]
        if not unflagged:
            self.log.warning("No unflagged antenna found to use as reference.")
            return None

        pos_array = np.array([pos for _, pos in unflagged])
        centroid  = pos_array.mean(axis=0)
        distances = np.linalg.norm(pos_array - centroid, axis=1)
        best      = unflagged[np.argmax(distances)][0]
        self.log.debug(f"Selected reference antenna (most distant from centre): {best} "
                      f"({distances.max():.0f} m from centroid)")
        return best

    def swap_feeds(self, column='DATA'):
        """
        swaps feeds of column (column) for the whole MS to correct for the wrong interpretation of the polarisation.
        This is needed for MeerKAT to have the correct polarisation information in the final products.

        Args:
            fields (list): list of field IDs to correct

        Returns:

        """        

        # Change RECEPTOR_ANGLE : DEFAULT IS -90DEG but should be fixed with the initial swap and set to 0 for the correct interpretation of the polarisation.
        with open_table(self.msfile+'::FEED') as tb:
            feed_angle = tb.getcol('RECEPTOR_ANGLE')
            if np.all(feed_angle == 0):
                self.log.info("Feeds already flipped, skipping.")
                return

        self.log.info("Flipping feeds for datacolumn: '{}'".format(column))

        with open_table(self.msfile+"::DATA_DESCRIPTION") as t:
            spwsel = t.getcol("SPECTRAL_WINDOW_ID")[0]
            poldescsel = t.getcol("POLARIZATION_ID")[0]

        with open_table(self.msfile+"::POLARIZATION") as t:
                poltype = t.getcol("CORR_TYPE").T[poldescsel] # trasnpose is only necessary if using CASA tables, not if using pyrap tables
                # must be linear
                if any(poltype - np.array([9,10,11,12]) != 0):
                    raise RuntimeError("Must be full correlation linear system being corrected")

        with open_table(self.msfile+"::SPECTRAL_WINDOW") as t:
                chan_freqs = t.getcol("CHAN_FREQ").T[spwsel]
                nchan = chan_freqs.size
            
            #timepaunix = np.array(list(map(lambda x: x.replace(tzinfo=pytz.UTC).timestamp(), timepadt)))
        nrowsput = 0
        with open_table(self.msfile, readonly=False) as t:
            nrow = t.nrows()
            nchunk = nrow // 1000 + int(nrow % 1000 > 0)
            for ci in range(nchunk):
                cl = ci * 1000
                crow = min(nrow - ci * 1000, 1000)
                data = t.getcol(column, startrow=cl, nrow=crow)  # shape: (4, nchan, crow)
                if data.shape[0] != 4:
                    raise RuntimeError("Data must be full correlation")
                data = data.T.reshape(crow, nchan, 2, 2)  # -> (crow, nchan, 2, 2)

                def give_crossphase_mat(phase, nrow, nchan, conjugate=False):
                    ones = np.ones(nchan*nrow)
                    zeros = np.zeros(nchan*nrow)
                    e = np.exp((1.j if not conjugate else -1.j) * np.deg2rad(phase)) * ones
                    return np.array([e,zeros,zeros,ones]).T.reshape(nrow, nchan, 2, 2)

                # need to apply anti-diagonal
                FVmat = np.array([np.zeros(nchan*crow),
                                            np.ones(nchan*crow),
                                            np.ones(nchan*crow),
                                            np.zeros(nchan*crow)]).T.reshape(crow, nchan, 2, 2)
                        
                # cojugate exp for left antenna
                XA1 = give_crossphase_mat(0.0, nrow=crow, nchan=nchan,conjugate=True)
                XA2 = give_crossphase_mat(0.0, nrow=crow, nchan=nchan,conjugate=False)

                JA1 = np.matmul(FVmat, XA1)
                JA2 = np.matmul(XA2, FVmat)

                corr_data = np.matmul(JA1, np.matmul(data, JA2)).reshape(crow, nchan, 4).T  # -> (4, nchan, crow)
                t.putcol(column, corr_data, startrow=cl, nrow=crow)
                self.log.debug("\tCorrected chunk {}/{}".format(ci+1, nchunk))
                nrowsput += crow
            assert nrow == nrowsput

        with open_table(self.msfile+'/FEED', readonly=False) as tb:
            feed_angle = tb.getcol('RECEPTOR_ANGLE')
            new_feed_angle = np.zeros(feed_angle.shape)
            tb.putcol('RECEPTOR_ANGLE', new_feed_angle)
