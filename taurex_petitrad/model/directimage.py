from .petitradtrans import petitRADTRANSModel
import numpy as np
from taurex.exceptions import InvalidModelException
from taurex.core import fitparam

class DirectImageRADTRANS(petitRADTRANSModel):

    @classmethod
    def input_keywords(self):
        return ['directimage-petitrad', 'direct-petitrad', ]

    

    def build_atmosphere_object(self):

        return self._radtrans.Radtrans(line_species=self.linespecies, \
                    rayleigh_species = self.rayleigh_species, \
                    continuum_opacities = self.continuum_species, \
                    wlen_bords_micron = self._wlen_micron)

    def path_integral(self, wngrid, return_contrib):
        from taurex.constants import RJUP,RSOL
        import astropy.units as u

        abundances, temperature, MMW, Rp, gravity, p0bar = self.setup_parameters()


        self.info('Molecular abundances at surface: %s',[ (k,v[-1]) for k,v in abundances.items()])
        self.info('Temperature at surface %s',temperature[-1])
        self.info('MMw at surface %s',MMW[-1])
        self.info('Planet radius: %s',Rp)
        self.info('Gravity in cm/2 at surface: %s',gravity)
        self.info('P0 = radius: %s',p0bar)



        self._atmosphere.calc_flux(temperature, abundances, gravity, MMW)
        integral = self._atmosphere.flux[::-1]
        petit_flux_u = integral * u.erg / u.cm**2 / u.s / u.Hz

        petit_flux_W = petit_flux_u.to(u.W/u.m**2/u.um, equivalencies=u.spectral_density(self.nativeWavenumberGrid*u.k)).value

        #print(integral)
        native = self.nativeWavenumberGrid
        native_filt = (native >= wngrid.min()) & (native <= wngrid.max())

        petit_flux_W = petit_flux_W[native_filt]
        if np.any(np.isnan(petit_flux_W)):
            raise InvalidModelException

        return petit_flux_W, np.zeros(shape=(self.nLayers,wngrid.shape[0]))