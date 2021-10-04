from .petitradtrans import petitRADTRANSModel
import numpy as np
from taurex.exceptions import InvalidModelException
from taurex.core import fitparam

class TransmissionRADTRANS(petitRADTRANSModel):

    @classmethod
    def input_keywords(self):
        return ['transmission-petitrad', 'transit-petitrad', ]

    def __init__(self, petitrad_path=None,
                    planet=None,
                    star=None,
                    pressure_profile=None,
                    temperature_profile=None,
                    chemistry=None,
                    nlayers=100,
                    atm_min_pressure=1e-4,
                    atm_max_pressure=1e6,
                    rayleigh_species = [],
                    continuum_species = [],
                    wlen_bords_micron=[0.3,15],
                    opacity_method=None,
                    Pcloud=None,
                    gamma_scat=None,
                    kappa_zero=None,
                    haze_factor = None,
                    include_condensates=False):
        super().__init__(petitrad_path,
                 planet,
                 star,
                 pressure_profile,
                 temperature_profile,
                 chemistry,
                 nlayers,
                 atm_min_pressure,
                 atm_max_pressure,
                 rayleigh_species,
                 continuum_species,
                 wlen_bords_micron,
                 opacity_method)

        self._Pcloud=None
        self._gamma_scat=gamma_scat
        self._kappa_zero=kappa_zero
        self._haze_factor = haze_factor
        self._cloud_pressure = Pcloud
        self.include_condensates = False

    def build_atmosphere_object(self):
        cloud_species = None
        if self.include_condensates:
            pass

        return self._radtrans.Radtrans(line_species=self.linespecies, \
                    rayleigh_species = self.rayleigh_species, \
                    continuum_opacities = self.continuum_species, \
                    wlen_bords_micron = self._wlen_micron)
    

    def path_integral(self, wngrid, return_contrib):
        from taurex.constants import RJUP,RSOL
        from petitRADTRANS import nat_cst as nc
        import astropy.units as u
        abundances, temperature, MMW, Rp, gravity, p0bar = self.setup_parameters()


        self.info('Molecular abundances at surface: %s',[ (k,v[-1]) for k,v in abundances.items()])
        self.info('Temperature at surface %s',temperature[-1])
        self.info('MMw at surface %s',MMW[-1])
        self.info('Planet radius: %s',Rp)
        self.info('Gravity in cm/2 at surface: %s',gravity)
        self.info('P0 = radius: %s',p0bar)

        Pcloud = self._cloud_pressure

        if Pcloud is not None:
            Pcloud = Pcloud*1e-5


        self._atmosphere.calc_transm(temperature, abundances, gravity, MMW, R_pl=Rp, P0_bar=p0bar, Pcloud=Pcloud,
                                     gamma_scat=self._gamma_scat,kappa_zero=self._kappa_zero, haze_factor=self._haze_factor,
                                     variable_gravity=True)

        Rs = self.star.radius*100
        integral = self._atmosphere.transm_rad**2
        # petit_flux_u = integral*u.erg/u.cm**2/u.s/u.Hz

        # petit_flux_W = petit_flux_u.to(u.W/u.m**2/u.um, equivalencies=u.spectral_density(self.nativeWavenumberGrid*u.k)).value

        #print(integral)

        rprs2 = (integral[::-1])/Rs**2

        native = self.nativeWavenumberGrid

        native_filt = (native >= wngrid.min()) & (native <= wngrid.max())

        rprs2 = rprs2[native_filt]
        if np.any(np.isnan(rprs2)):
            raise InvalidModelException

        return rprs2, np.zeros(shape=(self.nLayers,wngrid.shape[0]))

    @fitparam(param_name='kappa_zero', param_latex='$\kappa_0$',default_fit=False,default_mode='linear',default_bounds=[0.01,2.0])
    def kappaZero(self):
        return self._kappa_zero
    
    @kappaZero.setter
    def kappaZero(self, value):
        self._kappa_zero = value

    @fitparam(param_name='gamma_scat', param_latex='$\gamma$',default_fit=False,default_mode='linear',default_bounds=[-4,2])
    def gamma(self):
        return self._gamma
    
    @gamma.setter
    def gamma(self, value):
        self._gamma = value

    @fitparam(param_name='clouds_pressure',
              param_latex='$P_\mathrm{clouds}$',
              default_mode='log',
              default_fit=False, default_bounds=[1e-3, 1e6])
    def cloudsPressure(self):
        """
        Cloud top pressure in Pascal
        """
        return self._cloud_pressure

    @cloudsPressure.setter
    def cloudsPressure(self, value):
        self._cloud_pressure = value

    @fitparam(param_name='haze_factor',
              param_latex='$f_\mathrm{haze}$',
              default_mode='linear',
              default_fit=False, default_bounds=[1e-3, 1e6])
    def hazeFactor(self):
        """
        Cloud top pressure in Pascal
        """
        return self._haze_factor

    @hazeFactor.setter
    def hazeFactor(self, value):
        self._haze_factor = value

    #KZ11 KZ33 MK29

