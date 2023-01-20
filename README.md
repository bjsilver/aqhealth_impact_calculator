# Air quality health impact calculator

These scripts are used to calculate the health impact of air pollution
using
* modelled air pollution concentrations
* gridded population data
* population age structure data
* baseline health data


## Current limitations

* The scripts only use the GEMM health functions from Burnett et al. (2018)
but could later be expanded to use other health functions and other pollutants (e.g. O<sub>3</sub>)
* Currently setup to only use GPWv4 gridded population data
* Currently uses GBD IHME baseline health data and age structure

## Preparing to run the code

Download the gridded population data, baseline health data and age structure data

### Make population count slices
For the GPWv4 data, a utility script (`./utils/make_popslices.py`) has been written.
If using other gridded population data, you will need to make these files/edit the script.

For each year you wish to include in the analysis, there needs to be a pickled
dictionary of population count slices in the grids directory 
(e.g. `./grids/pop_count_slices_YYYY.P`)

The dictionary keys should be the 3 letter isocode for each country,
and the value should be a DataArray of population count for that country
with surrounding grids as zeroes e.g.

```
{'CHN': <xarray.DataArray 'CHN' (latitude: 860, longitude: 1469)>
 array([[0., 0., 0., ..., 0., 0., 0.],
        [0., 0., 0., ..., 0., 0., 0.],
        [0., 0., 0., ..., 0., 0., 0.],
        ...,
        [0., 0., 0., ..., 0., 0., 0.],
        [0., 0., 0., ..., 0., 0., 0.],
        [0., 0., 0., ..., 0., 0., 0.]])
 Coordinates:
   * longitude  (longitude) float64 73.6 73.65 73.69 73.73 ... 134.7 134.7 134.8
   * latitude   (latitude) float64 53.52 53.48 53.44 53.4 ... 16.06 16.02 15.77
 Attributes:
     units:      Persons
     long_name:  Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 a...
     min:        [ 0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.00...
     max:        [9.63072875e+05 1.00347762e+06 1.13291988e+06 1.34768400e+06\...}
```

### Make age structure lookup
For the IHME age structure data a utility script (`./utils/make_age_structures.py`) has been written.
If using other age structyre data, you will need to make these files/edit the script.

For the `PopulationData` class to be able to read in the age structure, it should
be in the format of a pickled pandas dataframe saved at `./lookups/age_structure_YYYY.P`

The dataframe should have a 3 level MultiIndex of 3 letter country isocode, age group,
and uncertainty level, with each value the proportion that the age group makes up in
the country's population, e.g.

```
country_isocode  age_group_name  uncertainty
ARM              25-29           val            0.078126
                                 upper          0.084380
                                 lower          0.072182
                 30-34           val            0.068297
                                 upper          0.073765
  
STP              75-79           upper          0.009699
                                 lower          0.008232
                 80+             val            0.008063
                                 upper          0.008764
                                 lower          0.007438
Length: 7344, dtype: float64
```

### Baseline health data
The code is currently written to run with the IHME baseline health data.
In the `hia_calculation.py` script, the class `BaseLineHealthData` is used to
read in the raw IHME csv as from the path defined in the config file parameter `bh_fpath`.

If using other baseline health data, it should be transformed to match the format
of the IHME data OR the `BaseLineHealthData` class should be rewritten.

## Running the code
First create a copy of the `config_template.yaml` and replace the settings/paths
with your own. Then the HIA should be submitted by running `hia.py`[1].

The code will first create a common grid based on the gridded population data
Then it will calculate the health impact for each country, age group, cause, and uncertainty level.

The results will be saved in a directory named according to the `scenario_name`
parameter in the config file. Each results directory should contain:
* by_country_results.csv
* by_country_age-group_uncertainty_results.csv  
* gridded_results.nc
* [scnario_name].yml (a copy of the config file)

[^1]: I doubt it's going to work though good luck