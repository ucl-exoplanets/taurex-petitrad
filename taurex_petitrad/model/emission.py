from .directimage import DirectImageRADTRANS

class EmissionRADTRANS(DirectImageRADTRANS):

    @classmethod
    def input_keywords(self):
        return ['emission-petitrad', 'eclipse-petitrad', ]
    

    def compute_final_flux(self, f_total):
        star_sed = self._star.spectralEmissionDensity

        self.debug('Star SED: %s', star_sed)
        # quit()
        star_radius = self._star.radius
        planet_radius = self._planet.fullRadius
        self.debug('star_radius %s', self._star.radius)
        self.debug('planet_radius %s', self._star.radius)
        last_flux = (f_total/star_sed) * (planet_radius/star_radius)**2

        self.debug('last_flux %s', last_flux)

        return last_flux

    def path_integral(self, wngrid,return_contrib):

        flux,tau = super().path_integral(wngrid,return_contrib)

        return self.compute_final_flux(flux),tau


