# HTCondor Lee-Wick leptons signal generation

This is a package that contains all the necessary ingredients to run the full `CMSSW` simulation of a Lee-Wick signal containing Lee-Wick (LW) electrons.  Additional theoretical details can be found elsewhere.

## Simulation data

Simulations satrt out from the *lhe* files stored in the `data` directory.  These were obtained beforehand using Feynrules and Madgraph.  This `data` directory also contains a text file called `VidasMedias.txt`. In this file we store:

`process spin charge color mass mass-width minmass maxmass halflife`

Only four mass points are currently implemented: 200, 300, 400 and 500 GeV for the LW electrons. 


The other folder in this package is the `config` directory.  Here one can find templates for the full simulation chain in `CMSSW`, namely *gensim*, *hlt* and *reco* sequences.  These were obtained using the *cmsDriver* command as instructed in the `cmsDriverLines.txt`, which was obtained following standard recipes.

## Simulation procedure

The scripts in this package are designed to run for HTCondor using continerized application.  In this particular setup, *Singularity* is used.

The `submit_jobs.sh` is the steering file.  This script loops over all the possible datasets to be simulated, creates the *jobs* using the `create_job.py` script and submits them to the Condor cluster.  The python script, `create_job.py` takes care of configuring the `job.sh` correctly so it can trigger a correct template for each mass point simulation.

Certain parameters in the headers of these scripts can be changed according to needs.
 

