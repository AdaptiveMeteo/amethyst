#!/usr/bin/env python

"""
RTTOV forward model adapter.

Implements the same compute(indata, outdata) protocol as ossFM so that
ForwardModel.py can use either OSS or RTTOV transparently.

indata convention (set by ForwardModel.py, all profiles top-first):
    temp     [nlev]        temperature, K
    h2o      [nlev]        water vapour, kg/kg
    o3       [nlev]        ozone, kg/kg
    co2      [1 or nlev]   CO2, kg/kg
    pressure [nlev]        pressure, hPa
    tskin    scalar        skin temperature, K
    psf      scalar        surface pressure, hPa
    sfgrd    [nsfgrd]      emissivity grid wavenumbers, cm-1
    emrf     [nsfgrd, 2]   emissivity (column 0 used)
    obsang   scalar        satellite zenith angle, degrees
    solzenith scalar       solar zenith angle, degrees
    azangle  scalar        solar azimuth angle, degrees
    lat      scalar        latitude, degrees

outdata keys filled (same layout as ossFM):
    y        [nchan]            radiances
    xkt      [4*nlev+2, nchan]  Jacobians (T | SKT | psf | WV | CO2 | O3)
                                    All from K-matrix for LM consistency.
    paxkemrf [2, nchan]         surface emissivity Jacobians (row 0 used)
    xkemrf   [1, nchan, 2]      placeholder (unused by ForwardModel.py)
"""

import numpy as np


class rttovFM(object):
    """RTTOV K-matrix wrapper with the ossFM compute() interface."""

    def __init__(self, coef_file, rttov_path, surftype=1, nthreads=1,
                 lib_path=None):
        """
        Parameters
        ----------
        coef_file : str
            Path to RTTOV coefficient file (.dat or .nc).
        rttov_path : str
            RTTOV installation root; wrapper/ subdirectory is added to sys.path.
        surftype : int
            RTTOV surface type: 0=land, 1=sea, 2=sea ice (default 1).
        nthreads : int
            RTTOV internal threads; keep 1 when using Python multiprocessing.
        lib_path : str, optional
            Extra directory prepended to LD_LIBRARY_PATH before importing
            pyrttov (needed when RTTOV was compiled in a different conda env).
        """
        import os, sys
        if lib_path:
            os.environ['LD_LIBRARY_PATH'] = (
                lib_path + ':' + os.environ.get('LD_LIBRARY_PATH', '')
            )
        sys.path.insert(0, rttov_path + '/wrapper')
        import pyrttov
        self._pyrttov = pyrttov

        rttov = pyrttov.Rttov()
        rttov.FileCoef = coef_file
        rttov.Options.CO2Data        = True
        rttov.Options.O3Data         = True
        rttov.Options.StoreRad       = True
        rttov.Options.ADKBT          = False   # K-matrix in radiance units (not BT)
        rttov.Options.VerboseWrapper = False
        rttov.Options.Verbose        = False   # suppress coefficient-limit warnings
        rttov.Options.ApplyRegLimits = True    # clamp profiles to valid coef range
        rttov.Options.CheckProfiles  = False   # disable hard-stop on unphysical values
        rttov.Options.Nthreads       = nthreads
        rttov.loadInst()   # load all channels

        self._rttov     = rttov
        self.cwvn       = np.array(rttov.WaveNumbers, dtype=np.float64)
        self.nchan      = len(self.cwvn)
        self._surftype  = surftype
        self._nsurfaces = 1

    # ------------------------------------------------------------------
    def compute(self, indata, outdata, dbg=False):
        """Run RTTOV K-matrix for one profile and fill outdata."""
        pyrttov = self._pyrttov
        rttov   = self._rttov

        temp = np.asarray(indata['temp'], dtype=np.float64)
        nlev = len(temp)
        nprof = 1
        nsurf = self._nsurfaces

        # --- Profiles -------------------------------------------------
        # Use ppmv_dry throughout so we can clip to RTTOV's documented
        # physical limits (max H2O ~500000 ppmv, max O3 ~50 ppmv).
        # During LM iterations xhat can push kg/kg values to astronomically
        # large numbers; RTTOV's internal physical-check aborts on those.
        Md = 28.966; Mw = 18.016; Mc = 44.01; Mo = 48.0

        myProf = pyrttov.Profiles(nprof, nlev + 1, nsurf)
        myProf.GasUnits = pyrttov.gasUnitType('ppmv_dry')

        myProf.P = indata['pressure'].reshape(1, nlev).astype(np.float64)
        temp_clip = np.clip(np.nan_to_num(temp, nan=200.0, posinf=300.0, neginf=100.0),
                            90.0, 400.0)   # also guards against large finite xhat values
        myProf.T = temp_clip.reshape(1, nlev)

        # H2O: kg/kg → ppmv_dry, clip to RTTOV physical bound.
        # Keep the raw (nan-sanitised but otherwise unclipped) value for the
        # Jacobian clip-ratio: wv_scale = q_clip/q_raw drives the effective
        # K toward zero when q diverges, giving correct chain-rule behaviour
        # in ForwardModel's log-space multiply (K_eff * vmr → finite).
        # Clip to 1 kg/kg ONLY for the RTTOV forward-model input to prevent
        # overflow in the ppmv multiply for large-but-finite xhat values.
        h2o_kgkg_raw = np.nan_to_num(np.asarray(indata['h2o'], dtype=np.float64),
                                      nan=1e-9, posinf=1e-2, neginf=1e-9)
        h2o_kgkg = np.clip(h2o_kgkg_raw, 0.0, 1.0)   # RTTOV input only
        h2o_ppmv = np.clip(h2o_kgkg * (Md / Mw) * 1e6, 1e-3, 500000.0)
        myProf.Q = h2o_ppmv.reshape(1, nlev)

        # O3: kg/kg → ppmv_dry, clip to physical bound (same strategy as H2O)
        o3_kgkg_raw = np.nan_to_num(np.asarray(indata['o3'], dtype=np.float64),
                                     nan=1e-9, posinf=1e-5, neginf=1e-9)
        o3_kgkg  = np.clip(o3_kgkg_raw, 0.0, 1e-2)   # RTTOV input only
        o3_ppmv  = np.clip(o3_kgkg  * (Md / Mo) * 1e6, 1e-6, 50.0)
        myProf.O3 = o3_ppmv.reshape(1, nlev)

        # CO2: already kg/kg from ForwardModel; convert to ppmv_dry
        co2 = np.asarray(indata['co2'], dtype=np.float64).ravel()
        if co2.size == 1:
            co2 = np.full(nlev, co2[0])
        co2_ppmv = np.clip(co2 * (Md / Mc) * 1e6, 100.0, 1000.0)
        myProf.CO2 = co2_ppmv.reshape(1, nlev)

        myProf.PHalf = self._build_phalf(
            indata['pressure'],
            float(np.asarray(indata['psf']).ravel()[0])
        ).reshape(1, nlev + 1)

        # --- Geometry -------------------------------------------------
        satzen = float(np.asarray(indata['obsang']).ravel()[0])
        sunzen = float(np.asarray(indata['solzenith']).ravel()[0])
        sunazi = float(np.asarray(indata['azangle']).ravel()[0])
        lat    = float(np.asarray(indata['lat']).ravel()[0])

        myProf.Angles    = np.array([[satzen, 0.0, sunzen, sunazi]], dtype=np.float64)
        myProf.SurfGeom  = np.array([[lat,    0.0, 0.0]],           dtype=np.float64)
        myProf.DateTimes = np.array([[2020, 1, 1, 0, 0, 0]],        dtype=np.int32)

        # --- Surface --------------------------------------------------
        tskin = float(np.asarray(indata['tskin']).ravel()[0])
        # Skin: [skin_T, salinity, snow_frac, foam_frac, fastem_coef x5]
        myProf.Skin     = np.array([[[tskin, 35., 0., 0., 3., 5., 15., 0.1, 0.3]]],
                                   dtype=np.float64)
        myProf.SurfType = np.array([[[self._surftype, 0]]], dtype=np.int32)

        # 2-metre variables: use lowest model level as proxy.
        # NearSurface q2m must match the gas unit (ppmv_dry).
        myProf.NearSurface = np.array(
            [[[float(temp_clip[-1]),
               float(h2o_ppmv[-1]),
               0., 0., 100000.]]],
            dtype=np.float64
        )

        rttov.Profiles = myProf

        # --- Emissivity -----------------------------------------------
        # Interpolate AMETHYST eigenvalue-based emissivity to RTTOV channels.
        # Positive values tell RTTOV to use the supplied emissivity directly.
        sfgrd = np.asarray(indata['sfgrd']).ravel()
        emis  = np.asarray(indata['emrf'])[:, 0]
        emis_rttov = np.interp(self.cwvn, sfgrd, emis,
                               left=emis[0], right=emis[-1])
        surfemisrefl = np.zeros((5, nprof, nsurf, self.nchan), dtype=np.float64)
        surfemisrefl[0, 0, 0, :] = emis_rttov
        rttov.SurfEmisRefl = surfemisrefl

        # --- Run K-matrix ---------------------------------------------
        rttov.runK()

        # Extract all K-matrix results immediately.
        y0   = np.array(rttov.Rads[0, :], dtype=np.float64)
        TK   = np.array(rttov.TK[0],       dtype=np.float64)
        QK   = rttov.getItemK('Q')
        CO2K = rttov.getItemK('CO2')
        O3K  = rttov.getItemK('O3')
        _SKK = rttov.SkinK     # shape (nprof, nsurf, nchan, 9); index 0 = T component
        _SEK = rttov.SurfEmisK # shape (nprof, nsurf, nchan) or None

        # --- Radiances ------------------------------------------------
        outdata['y'] = y0

        # --- Jacobians in OSS xkt layout ------------------------------
        # Row layout: 0:nlev=T, nlev=SKT, nlev+1=psf(unused),
        #             nlev+2:2n+2=WV, 2n+2:3n+2=CO2, 3n+2:4n+2=O3
        xkt = np.zeros((4 * nlev + 2, self.nchan), dtype=np.float64)

        # T: TK[0] shape (nchan, nlev) in pyrttov
        xkt[0:nlev, :] = TK.T if TK.shape == (self.nchan, nlev) else TK

        # SKT: SkinK index 0 = temperature component; in radiance units (ADKBT=False)
        if _SKK is not None:
            xkt[nlev, :] = np.array(_SKK[0, 0, :, 0], dtype=np.float64)

        # WV — RTTOV gives dF/d(q_ppmv_dry); ForwardModel.py computes
        # dF/d(log q) = K_Q_kgkg * exp(xhat_wv).  We must account for
        # any clipping: K_Q_eff = K_Q_ppmv * (Md/Mw)*1e6 * (q_clip/q_raw).
        # Use h2o_kgkg_raw (= vmr from ForwardModel, before the RTTOV clip)
        # so that wv_scale = q_clip/q_raw → 0 when vmr diverges, which
        # cancels ForwardModel's *vmr multiply and keeps K finite.
        q_clip_kgkg = h2o_ppmv * (Mw / Md) * 1e-6   # kg/kg actually used by RTTOV
        wv_scale = q_clip_kgkg / np.maximum(h2o_kgkg_raw, 1e-30)  # use raw value
        if QK is not None:
            QK0 = np.array(QK[0], dtype=np.float64)
            QK0 = QK0.T if QK0.shape == (self.nchan, nlev) else QK0
            xkt[nlev + 2:2 * nlev + 2, :] = QK0 * (Md / Mw) * 1e6 * wv_scale[:, None]

        # CO2 — ppmv→kg/kg conversion; CO2 is typically unclipped so scale ≈ 1
        if CO2K is not None:
            CO2K0 = np.array(CO2K[0], dtype=np.float64)
            CO2K0 = CO2K0.T if CO2K0.shape == (self.nchan, nlev) else CO2K0
            xkt[2 * nlev + 2:3 * nlev + 2, :] = CO2K0 * (Md / Mc) * 1e6

        # O3 — same clip-ratio correction as WV (use raw value for same reason)
        o3_clip_kgkg = o3_ppmv * (Mo / Md) * 1e-6
        o3_scale = o3_clip_kgkg / np.maximum(o3_kgkg_raw, 1e-30)  # use raw value
        if O3K is not None:
            O3K0 = np.array(O3K[0], dtype=np.float64)
            O3K0 = O3K0.T if O3K0.shape == (self.nchan, nlev) else O3K0
            xkt[3 * nlev + 2:4 * nlev + 2, :] = O3K0 * (Md / Mo) * 1e6 * o3_scale[:, None]

        outdata['xkt'] = xkt

        # --- Surface emissivity Jacobians ----------------------------
        if 'paxkemrf' not in outdata or outdata['paxkemrf'].shape[1] != self.nchan:
            outdata['paxkemrf'] = np.zeros((2, self.nchan), dtype=np.float64)
        if _SEK is not None:
            outdata['paxkemrf'][0, :] = np.array(_SEK[0, 0, :], dtype=np.float64)

        # Placeholder kept for API compatibility with ossFM
        if 'xkemrf' not in outdata:
            outdata['xkemrf'] = np.zeros((1, self.nchan, 2), dtype=np.float64)

    # ------------------------------------------------------------------
    @staticmethod
    def _build_phalf(p_full, p_surface):
        """Half-level pressures from full levels and surface pressure.

        RTTOV requires strict interleaving: ph[i] < p[i] < ph[i+1].
        The bottom half-level must be strictly below the last full level,
        which is not guaranteed when AMETHYST's profile touches the surface.
        """
        p = np.asarray(p_full, dtype=np.float64)
        nlev = len(p)
        ph = np.empty(nlev + 1, dtype=np.float64)
        ph[0]    = max(1.0e-6, 0.5 * p[0])
        ph[1:-1] = 0.5 * (p[:-1] + p[1:])
        # Extrapolate surface half-level one half-interval below the last
        # full level, or use p_surface — whichever is larger.
        dp_last = p[-1] - p[-2]
        ph[-1]   = max(p_surface, p[-1] + 0.5 * dp_last)
        return ph
