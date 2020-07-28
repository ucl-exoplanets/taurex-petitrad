from taurex.opacity import InterpolatingOpacity
import pathlib
import numpy as np
import os
import pathlib
from taurex.util.util import calculate_weight

class LineByLine(InterpolatingOpacity):

    @classmethod
    def priority(cls):
        return 200

    @classmethod
    def discover(cls):
        import os
        import glob
        import pathlib
        from taurex.cache import GlobalCache
        from taurex.util.util import sanitize_molecule_string

        input_path = None
        if GlobalCache()['petitrad_path'] is not None:
            input_path = os.path.join(GlobalCache()['petitrad_path'],'petitRADTRANS','input_data')

        input_path = input_path or GlobalCache()['xsec_path']
        if input_path is None:
            return []
        opacities_path = os.path.join(input_path,'opacities','lines','line_by_line')
        if not os.path.isdir(opacities_path):
            return []
        
        folders = [x for x in glob.glob(os.path.join(opacities_path,'*')) if os.path.isdir(x)] 

        check_molecules = zip([pathlib.Path(x).stem.split('_',1) for x in folders],folders)

        molecule_dict = {}

        main_iso_molecules = [(a,b) for a,b in check_molecules if a[-1]=='main_iso']

        for mol,filename in main_iso_molecules:
            molecule_dict[mol[0]] = filename

        check_molecules = zip([pathlib.Path(x).stem.split('_',1) for x in folders],folders)

        for mol,filename in check_molecules:
            if len(mol) == 1 and mol[0] not in molecule_dict:
                molecule_dict[mol[0]] = filename 


        return list(molecule_dict.items())


    def __init__(self, molecule_path, interpolation_mode='linear'):
        super().__init__('[PR]LineByLineOpacity:{}'.format(pathlib.Path(molecule_path).stem[0:10]),
                        interpolation_mode=interpolation_mode)

        self.info('Loading molecule from path %s', molecule_path)
        self._molecule_path = molecule_path
        self._molecule_name = pathlib.Path(molecule_path).stem.split('_',1)[0]
        
        self._determine_grids()
        self._load_xsec_from_path()

    def _load_xsec_from_path(self):
        self._wavenumber_grid = 10000/(np.fromfile(os.path.join(self._molecule_path,'wlen.dat'))[::-1]*1e4)
        self._xsec_grid = np.empty(shape=(self._pressure_grid.shape[0], 
                                          self._temperature_grid.shape[0],
                                          self._wavenumber_grid.shape[0]))

        num_moles = 1/calculate_weight(self.moleculeName)

        num_molecules = num_moles*6.0221409e23

        for p, t, sigma in self._sigma_files:
            pindex = np.where(self._pressure_grid==p)[0]
            tindex = np.where(self._temperature_grid==t)[0]



            self._xsec_grid[pindex, tindex,:] = np.fromfile(sigma)[::-1]/num_molecules
    
    def _determine_grids(self):
        import os
        import glob

        sigma_files = glob.glob(os.path.join(self._molecule_path,'sigma_*.dat'))
        temperature_grid = [float(pathlib.Path(x).stem.split('_')[2][:-1]) for x in sigma_files]
        pressure_grid = [float(pathlib.Path(x).stem.split('_')[3].split('b')[0])*1e5 for x in sigma_files]

        #print([pathlib.Path(x).stem.split('_')[3].split('b')[0] for x in sigma_files])

        s = sorted(zip(pressure_grid,temperature_grid,sigma_files),key = lambda x:(x[0],x[1]))
        pressure_grid,temperature_grid,sigma_files = list(zip(*s))

        self._pressure_grid = np.unique(np.array(pressure_grid))
        self._temperature_grid = np.unique(np.array(temperature_grid))
        self._sigma_files = list(zip(pressure_grid,temperature_grid,sigma_files))
        self._max_pressure = self._pressure_grid.max()
        self._min_pressure = self._pressure_grid.min()
        self._max_temperature = self._temperature_grid.max()
        self._min_temperature = self._temperature_grid.min()
    @property
    def wavenumberGrid(self):
        return self._wavenumber_grid

    @property
    def temperatureGrid(self):
        return self._temperature_grid

    @property
    def pressureGrid(self):
        return self._pressure_grid

    @property
    def moleculeName(self):
        return self._molecule_name

    @property
    def xsecGrid(self):
        return self._xsec_grid