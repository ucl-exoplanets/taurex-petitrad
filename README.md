# TauREx-petitRADTRANS plugin

Current version: 0.0.0-alpha

TauREx-petitRADTRANS is a plugin that allows for the use of [petitRADTRANS](https://petitradtrans.readthedocs.io/en/latest/) opacities
and forward models in [TauREx](https://github.com/ucl-exoplanets/TauREx3_public)

## Installing

Installing the plugin from PyPI requires only a single command:
```bash
pip install taurex_petitrad
```
You can also install from source by doing:
```bash
git clone https://github.com/ucl-exoplanets/taurex-petitrad
cd taurex-petitrad
pip install .
```

## Using opacities

Once installed you have access to both petitRADTRANS line-by-line cross-sections and k-tables. You can
either set the *xsec_path* or *ktable_path* to the **input_data** folder in petitRADTRANS:
```
[Global]
xsec_path = /path/to/petitRADTRANS/petitRADTRANS/input_data/
ktable_path = /path/to/petitRADTRANS/petitRADTRANS/input_data/
```

Or by setting *petitrad_path* to the root petitRADTRANS folder

```
[Global]
petitrad_path = /path/to/petitRADTRANS/petitRADTRANS/
```

You can include both if you desire, *petitrad_path* will be prioritized.

## Using forward models

To use the forward models, the *petitrad_path* must be set to the root petitRADTRANS folder.
TauREx 3 will automatically include the package into Python so there is no need to modify system paths.

Once done, you will have access to three new forward models:

## Transmission

Transmission spectra can be accessed through the *model_type* keyword by setting it
as either **transmission-petitrad** or **transit-petitrad**
```
[Model]
model_type = transmission-petitrad
continuum_opacities = H2-H2, H2-He
wlen_bords_micron = 0.3,15.0
```

The available arguments are:

|Argument| Description| Type| Default | Required |
---------|------------|-----|---------|----------|
rayleigh_species | Rayleigh species to include | list of molecules| None | |
continuum_species | CIA species | list of molecule pairs| None| |
wlen_bords_micron | Wavelength range in um | array| 0.3, 15| |
Pcloud | Grey cloud deck pressure in Pa | float | None| |
gamma_scat | Powerlaw index | integer| None | |
kappa_zero | Powerlaw cloud opacity term in cm2/g | float | None | |
haze_factor | Rayleigh opacity scale term | float | None | |


For the more through explaination of cloud terms
please refer to [this](https://petitradtrans.readthedocs.io/en/latest/content/notebooks/clouds.html) documentation.

The available retrieval parameters are:

|Fitting Parameter| Description| 
Pcloud | Grey cloud deck pressure in Pa |
gamma_scat | Powerlaw index |
kappa_zero | Powerlaw cloud opacity term in cm2/g |
haze_factor | Rayleigh opacity scale term |

## Eclipse and Direct imaging

The eclipse spectra can be generated using the **emission-petitrad** or **eclipse-petitrad** keyword
```
[Model]
model_type = eclipse-petitrad
continuum_opacities = H2-H2, H2-He
wlen_bords_micron = 0.3,15.0
```
And direct imaging as through by **directimage-petitrad** or **direct-petitrad**:
```
[Model]
model_type = directimage-petitrad
continuum_opacities = H2-H2, H2-He
wlen_bords_micron = 0.3,15.0
```

Both have the same arguments available:
|Argument| Description| Type| Default | Required |
---------|------------|-----|---------|----------|
rayleigh_species | Rayleigh species to include | list of molecules| None | |
continuum_species | CIA species | list of molecule pairs| None| |
wlen_bords_micron | Wavelength range in um | array| 0.3, 15| |

There are no extra retrieval parameters provided by these forward models

## Limitations

Currently no condensate opacities are supported as the current release of 
TauREx does not have a comprehensive enough scattering framework to prevent incompatibility with other
models. Future releases will allow for their inclusion.

The majority of TauREx and its plugins can be used with the petitRADTRANS models such as chemistries,
temeperature profiles, samplers, non-uniform priors etc.
The only exception are the *contributions*. This is because the integral is being performed
by petitRADTRANS itself, therefore some plugins (such as TauREx-CUDA) cannot be used. Attempting to add
any TauREx contribution object such as:
```
[Model]
model_type = transmission-petitrad
continuum_opacities = H2-H2, H2-He
wlen_bords_micron = 0.3,15.0
    [[Absorption]]

    [[HydrogenIon]]

    [[Rayleigh]]
```

Will result in the error:
```
ValueError: TauREx 3 contributions are not supported with the petitRADTRANS forward models
```


