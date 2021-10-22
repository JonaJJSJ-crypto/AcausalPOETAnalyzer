# HTCondor Lee-Wick leptons signal generation

This is a package that contains all the necessary ingredients to run the full `CMSSW` simulation of a Lee-Wick signal containing Lee-Wick (LW) electrons.  Additional theoretical details can be found elsewhere.

## Simulation data

Simulations stars out from the *lhe* files stored in the `data` directory. If they are tared and/or zipped, they need to be decompressed after cloning the repository. These were obtained beforehand using Feynrules and Madgraph.  This `data` directory also contains a text file called `VidasMedias.txt`. In this file we store:

`process spin charge color mass mass-width minmass maxmass halflife`

Only four mass points are currently implemented: 200, 300, 400 and 500 GeV for the LW electrons. 


The other folder in this package is the `config` directory.  Here one can find templates for the full simulation chain in `CMSSW`, namely *gensim*, *hlt* and *reco* sequences.  These were obtained using the *cmsDriver* command as instructed in the `cmsDriverLines.txt`, which was obtained following standard recipes.  There is also a configuration file that allows one to merge the AODSIM output files.  This is necessary if there are too many files (with too few events), which can make it difficult to process.

## Simulation procedure

The scripts in this package are designed to run for HTCondor using continerized application.  In this particular setup, *Singularity* is used.

The `submit_jobs.sh` is the steering file.  This script loops over all the possible datasets to be simulated, creates the *jobs* using the `create_job.py` script and submits them to the Condor cluster.  The python script, `create_job.py` takes care of configuring the `job.sh` correctly so it can trigger a correct template for each mass point simulation.

Certain parameters in the headers of these scripts can be changed according to needs.

Finally, a `checkFiles.py` script is used to check for the availability of the resulting ROOT files of the simulation.  If some files are missing a `failedarguments` arguments and jdl files are created so one can resend that specific job by hand. 

## Meging AODSIM

After all jobs have been successfully executed, one can merge these AODSIM files into a reduced number of files (with increased size).  This could help with the further processing of the data.  For this, similar scripts (designed for Condor submission) are available.  These contain the word `merge` in order to differentiate them.  There are some parameters that need to be adjusted in the headers of these files according to the needs.


