from .petitradtrans import petitRADTRANSModel
import numpy as np

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
                    wlen_bords_micron=[0.3,15]):
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
                 wlen_bords_micron)
        pass


    def build_atmosphere_object(self):
        return self._radtrans.Radtrans(line_species=self.linespecies, \
                    rayleigh_species = self.rayleigh_species, \
                    continuum_opacities = self.continuum_species, \
                    wlen_bords_micron = self._wlen_micron)
    

    def path_integral(self, wngrid, return_contrib):
        
        abundances, temperature, MMW, Rp, gravity, p0bar = self.setup_parameters()

        self._atmosphere.calc_transm(temperature, abundances, gravity, MMW, R_pl=Rp, P0_bar=p0bar)

        F_lambda = self._atmosphere.transm_rad/Rp

        Rs = self.star.radius*100

        rprs2 = (F_lambda[::-1]*Rp/Rs)**2

        return rprs2/2, np.zeros(shape=(self.nLayers,rprs2.shape[0]))






    #KZ11 KZ33 MK29

