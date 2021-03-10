
import numpy as np
from taurex.model import SimpleForwardModel
import sys, os
from taurex.cache import GlobalCache
from ..opacities.linebyline import LineByLine
from ..opacities.ck import CKTable
from taurex.core import fitparam
from ..util import to_mass_frac
class petitRADTRANSModel(SimpleForwardModel):

    imported = False


    def __init__(self,
                 petitrad_path = None,
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
                 P0=-1,
                 ):
        super().__init__(self.__class__.__name__,
                         planet,
                         star,
                         pressure_profile,
                         temperature_profile,
                         chemistry,
                         nlayers,
                         atm_min_pressure,
                         atm_max_pressure)
        if not self.imported:
            attempt_load = False
            try:
                import petitRADTRANS
            except (ImportError, ModuleNotFoundError,):
                attempt_load = True
            
            if attempt_load:
                path = petitrad_path or GlobalCache()['petitrad_path']
                GlobalCache()['petitrad_path'] = path

                if path is None:
                    self.error('No path found in either [Global] petitrad_path or ')
                    raise ValueError('No path to petitRADTRANS set in either [Global] or initialization')
                
                sys.path.append(path)
                try:
                    import petitRADTRANS
                    self.imported = True
                except (ImportError, ModuleNotFoundError,):
                    self.error('Could not import petitRADTRANS from path %s', path)
                    raise ImportError('Could not load petitRADTRANS module')

        import petitRADTRANS

        opacity_method = opacity_method or  GlobalCache()['opacity_method']
        GlobalCache()['opacity_method'] = opacity_method
        
        if opacity_method in ('ktables',):
            opacity_method = 'c-k'
        else:
            opacity_method = 'lbl'
        
        self._opacity_method = opacity_method

        self._wlen_micron = wlen_bords_micron

        self._radtrans = petitRADTRANS

        self._initialized_petit = False
        self.linespecies = None
        self.rayleigh_species = rayleigh_species
        self.continuum_species = continuum_species
        self._P0 = P0

    def initialize_profiles(self):
        from taurex.util.util import conversion_factor
        super().initialize_profiles()

        if not self._initialized_petit:
            self.info('Initializing petitRADTRANS')
            self.linespecies = self.get_linespecies()
            self.info('Detected active line species %s',self.linespecies)
            self._atmosphere = self.build_atmosphere_object()
            Pabar = conversion_factor('Pa','bar')
            pressures = self.pressureProfile[::-1]*Pabar
            self.info('Pressure at surface %s',pressures[-1])
            self._atmosphere.setup_opa_structure(pressures)
            self._initialized_petit = True
    def build_atmosphere_object(self):
        raise NotImplementedError


    def get_linespecies(self):
        import pathlib
        from taurex.util.util import get_molecular_weight
        opacity_list = None
        if self._opacity_method == 'c-k':
            opacity_list = CKTable.discover()
        elif self._opacity_method == 'lbl':
            opacity_list = LineByLine.discover()
        
        opacity_dictionary = {x:pathlib.Path(y).stem for x,y in opacity_list}
        active_gases = self.chemistry.activeGases

        return list(set([opacity_dictionary[x] for x in active_gases]))


    
    @property
    def nativeWavenumberGrid(self):
        from taurex.constants import SPDLIGT
        return self._atmosphere.freq[::-1]/SPDLIGT/100

    def build_abundance(self):
        import pathlib
        from taurex.util.util import get_molecular_weight
        from taurex.constants import AMU

        active_gases = self.chemistry.activeGases
        active_profile =self.chemistry.activeGasMixProfile[:,::-1]
        inactive_gases = self.chemistry.inactiveGases
        inactive_profile = self.chemistry.inactiveGasMixProfile[:,::-1]
        mu= self.chemistry.muProfile[::-1]/AMU

        abundance = {}

        for mol,prof,petit_name in zip(active_gases,active_profile,self.linespecies):
            abundance[petit_name] =to_mass_frac(mol,prof,mu)
        
        # for mol in self.rayleigh_species:
        #     if mol in abundance:
        #         continue
        #     if mol in active_gases or mol in self.chemistry.inactiveGases:
        #         mol_weight = get_molecular_weight(mol)
        #         prof = self.chemistry.get_gas_mix_profile(mol)[::-1]
        #         abundance[mol] = to_mass_frac(mol,prof,mu)

        # for cia in self.continuum_species:

        #     for mol in cia.split('-'):
        #         if mol in abundance:
        #             continue
        #         if mol in active_gases or mol in self.chemistry.inactiveGases:
        #             mol_weight = get_molecular_weight(mol)
        #             prof = self.chemistry.get_gas_mix_profile(mol)[::-1]
        #             abundance[mol] = to_mass_frac(mol,prof,mu)
        for mol,prof in zip(inactive_gases,inactive_profile):
            abundance[mol] =to_mass_frac(mol,prof,mu)

        return abundance

    def setup_parameters(self):
        from taurex.util.util import conversion_factor
        from taurex.constants import AMU
        from petitRADTRANS import nat_cst as nc

        Pabar = conversion_factor('Pa','bar')
        
        abundances = self.build_abundance()

        temperature = self.temperatureProfile[::-1]

        MMW = self.chemistry.muProfile[::-1]/AMU

        Rp = self.planet.radius*nc.r_jup_mean

        P0 = self._P0
        if self._P0 == -1:
            P0 = 1e999

        P0 = min(self.pressureProfile[0],P0)
        P0 = max(self.pressureProfile[-1],P0)

        p0bar = P0*Pabar

        idx = (np.abs(self.pressureProfile - P0)).argmin()

        altitude = self.altitudeProfile[idx]

        gravity = self.planet.gravity_at_height(altitude)*100


        return abundances, temperature, MMW, Rp, gravity, p0bar

    @fitparam(param_name='P0_radius', param_latex='P$_0$',default_fit=False,default_mode='log',default_bounds=[2,-5])
    def P0(self):
        return self._P0
    
    @P0.setter
    def P0(self, value):
        self._P0 = value


    def add_contribution(self, contribution):
        raise ValueError('TauREx 3 contributions are not supported with the petitRADTRANS forward models')

    BIBTEX_ENTRIES = [
        """
        @article{Mollire_petitrad,
            title={petitRADTRANS},
            volume={627},
            ISSN={1432-0746},
            url={http://dx.doi.org/10.1051/0004-6361/201935470},
            DOI={10.1051/0004-6361/201935470},
            journal={Astronomy \& Astrophysics},
            publisher={EDP Sciences},
            author={Mollière, P. and Wardenier, J. P. and van Boekel, R. and Henning, Th. and Molaverdikhani, K. and Snellen, I. A. G.},
            year={2019},
            month={Jul},
            pages={A67}
        }""",

    ]





        
