
import numpy as np
from taurex.model import SimpleForwardModel
import sys, os
from taurex.cache import GlobalCache
from ..opacities.linebyline import LineByLine
from ..opacities.ck import CKTable

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
                 wlen_bords_micron=[0.3,15]
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

        opacity_method = GlobalCache()['opacity_method']

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
            self._atmosphere.setup_opa_structure(pressures)
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

        return [opacity_dictionary[x] for x in active_gases]


    
    @property
    def nativeWavenumberGrid(self):
        from taurex.constants import SPDLIGT
        return self._atmosphere.freq[::-1]/SPDLIGT/100

    def build_abundance(self):
        import pathlib
        from taurex.util.util import get_molecular_weight

        active_gases = self.chemistry.activeGases
        active_profile =self.chemistry.activeGasMixProfile

        muProfile = self.chemistry.muProfile

        abundance = {}

        for mol,prof,petit_name in zip(active_gases,active_profile,self.linespecies):
            mol_weight = get_molecular_weight(mol)
            abundance[petit_name] = mol_weight*prof[::-1]/muProfile[::-1]
        
        return abundance

    def setup_parameters(self):
        from taurex.util.util import conversion_factor
        from taurex.constants import AMU
        Pabar = conversion_factor('Pa','bar')
        
        abundances = self.build_abundance()

        temperature = self.temperatureProfile[::-1]

        MMW = self.chemistry.muProfile/AMU

        Rp = self.planet.fullRadius*100

        gravity = self.planet.gravity*10

        p0bar = self.pressureProfile[0]*Pabar

        return abundances, temperature, MMW, Rp, gravity, p0bar





        
