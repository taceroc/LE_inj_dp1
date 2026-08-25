# Inject LE simulation on DP1[1] images

using the LSST Sciennce Pipeline[2] version `lsst_stack = v29.1.1` on NERSC.

## How to create the DP1 images?
* `create_coadds/dp1_extract_visits.ipynb`: 
	1. For DP1 data select a tract and a band.
		- https://dp1.lsst.io/overview/observations.html use it for reference
		- for my particular case I selected the tract=4848 and band=r in the ECDFS field
	2. Select two MJDs visits windows, where each window has a similar number of visits.
		- I selected [60623-60633] and [60651-60655]
	The notebook will print the id for the individual visits in the format needed for the next step.
* `create_coadds/dp1_build_coadd.sh`: creates the coadds by calling some parts of the step3 of the pipeline. Here you will need to copy and paste the output of the  `dp1_extract_visits`.


## How to inject the LE simulations on the coadds?

### First, generate the LE image, and save it as a fits file 

The [LEIMG](https://github.com/taceroc/LEIMG) repo on the main branch, contains the code to generate the LE.
(this doesn't need the lsst_stack)
```
python LEIMG/main.py SimulateLEInfPlane -file_to_parameters <path-to-input-yml> -outdir <output-folder>
```

You need the `output-folder` and their format for the LE injections.

### Second, inject the LE on the coadds

You have to load first the `lsst_stack`
```
export STACKCVMFS=/cvmfs/sw.lsst.eu/almalinux-x86_64/lsst_distrib
export LSST_STACK_VERSION=v29.1.1

module load cpu

source $STACKCVMFS/$LSST_STACK_VERSION/loadLSST-ext.bash
setup -t v29_1_1 lsst_distrib

export DAF_BUTLER_REPOSITORY_INDEX=/global/cfs/cdirs/lsst/production/gen3/shared/data-repos.yaml
```

The scripts do the injection of the simulated LE into DP1 coadds, do source detection, perfom difference imaging, and creates and saves the triplets postage stamps as `.npy` and `png` format.

#### There are 4 different approaches to inject LE into the coadds
*Everything is still work in progress, `.py` scripts still need manual inputs and modification*

**CASE A: LE is only observed on the first epoch**

`python inject_diff_save_dp1_case_a.py -path_to_manifest <output-folder>/manifest.yml` 
* You need to define where the coadds are located

The injections will be saved as:

```
output-folder
├── injections
│   ├── asingle
│   │   ├── images
│   │   ├── numpy
│	└── a_injections_single_save_runs.csv
```

**CASE C: LE is observed at two different epochs, and the difference has the characteristic *black and white* for the observed LE**

`python inject_diff_save_dp1_case_c.py -path_to_manifest <output-folder>/manifest.yml` 
* You need to define where the coadds are located

The injections will be saved as:

```
output-folder
├── injections
	├── asingle
│   ├── cdouble
│   │   ├── images
│   │   ├── numpy
│	└── c_injections_double_onesys_save_runs.csv
```

**CASE B: Two LEs are observed at the first epoch; that is the same for both**

**CASE D: Two LEs are observed at two different epochs**



[1] NSF-DOE Vera C. Rubin Observatory (2025); Legacy Survey of Space and Time Data Preview 1, https://doi.org/10.71929/rubin/2570308

[2] Rubin Observatory Science Pipelines Developers (2025); The LSST Science Pipelines Software: Optical Survey Pipeline Reduction and Analysis Environment, https://doi.org/10.71929/rubin/2570545


