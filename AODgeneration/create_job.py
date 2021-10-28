#!/usr/bin/env python


import os
import sys

#Assuming each LW mass point has 150000 events, 
# we set the number of jobs to 15000, so there are 10 events per job to complete 150000
# we do not do 
numJobs = 7500
#This will be picked up automatically by job.sh: 
evtsPerJob = 20
#any number to start:
theFirstRun = 200000

jdl = """\
executable = ./{PROCESS}.sh
output = out/$(ProcId).$(ClusterID).out
error = err/$(ProcId).$(ClusterID).err
log = log/$(ProcId).$(ClusterID).log
max_retries = 3
RequestCpus = 1
+MaxRuntime = 3600
queue arguments from arguments.txt\
"""


def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)


def parse_arguments():
    if not len(sys.argv) == 3:
        raise Exception("./create_job.py PROCESS PATH_TO_JOBDIR")
    return {"process": sys.argv[1], "jobdir": sys.argv[2]}


def main(args):
    process = args["process"]
    print("Process: %s" % process)

    #parse half lifes file with information about mass 
    if os.path.isfile('data/VidasMedias.txt'):
        for line in open('data/VidasMedias.txt').readlines():
            if process in line:
                pythiaconf = line.split()
                thespin = pythiaconf[1]
                print ("thespin = "+thespin) 
                thecharge = pythiaconf[2]
                print ("thecharge = "+thecharge) 
                thecolor = pythiaconf[3]
                print ("thecolor = "+thecolor) 
                themass = pythiaconf[4]
                print ("themass = "+themass) 
                themwidth = pythiaconf[5]
                print ("themwidth = "+themwidth) 
                theminmass = pythiaconf[6]
                print ("theminmass = "+theminmass) 
                themaxmass = pythiaconf[7]
                print ("themaxmass = "+themaxmass) 
                thetau = pythiaconf[8]
                print ("thetau = "+thetau) 

    # Build argument list 
    arguments = []
    counter = 0
    for counter in range(numJobs):
        skipEvents = counter*evtsPerJob
        firstEvent = skipEvents+1
        firstRun = theFirstRun+counter
        arguments.append("%u %s %u %u %u %u %s %s" % (counter, process, evtsPerJob, skipEvents, firstEvent, firstRun, themass, thetau))
    #print(arguments)    
    print("Number of jobs: %u" % len(arguments))

    # Create jobdir and subdirectories
    jobdir = os.path.join(args["jobdir"], process)
    print("Jobdir: %s" % jobdir)
    mkdir(jobdir)
    mkdir(os.path.join(jobdir, "out"))
    mkdir(os.path.join(jobdir, "log"))
    mkdir(os.path.join(jobdir, "err"))

    # Write jdl file
    out = open(os.path.join(jobdir, "job.jdl"), "w")
    out.write(jdl.format(PROCESS=process))
    out.close()

    # Write argument list
    arglist = open(os.path.join(jobdir, "arguments.txt"), "w")
    for a in arguments:
        arglist.write(a+"\n")
    arglist.close()

    # Write job file
    jobfile = open("job.sh", "r").read()
    job = open(os.path.join(jobdir, "{PROCESS}.sh".format(PROCESS=process)), "w")
    job.write(jobfile)
    job.close()


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
