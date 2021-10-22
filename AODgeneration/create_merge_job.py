#!/usr/bin/env python

import os
import sys

#Ideally there should be less than 400 files for each
#index file that will be merged.  This is because each file
#seems to take about 6 seconds to process, so keeping jobs with less than 400
#keeps jobs under an hour
#ALSO IMPORTANT is to keep number of events under 6K per job, otherwise
#with the file size limit, more than one file gets created and the extra ones
#are not going to be copied. 6.5K events give approximately reaches the limit
nlinesSplit = 375

mergeDirName = "merged"

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
    print("Merging process: %s" % process)
    indexFile = "data/CMSPrivate_MonteCarlo_"+process+"_unmerged_file_index.txt"
    print("Splitting file: %s" % indexFile)
    splitBaseName = "CMSPrivate_MonteCarlo_"+process+"_split_"
    #split the index file
    splitline = "split -l "+str(nlinesSplit)+" "+indexFile+" "+splitBaseName+" --additional-suffix='.txt'" 
    os.system(splitline)

 # Build argument list                                                                                                                                                                                       
    print("Filelist:")
    arguments = []
    counter = 0
    processSplit = process+"_split_"
    for filename in os.listdir("."):
        if processSplit in filename:
            print("    %s." % filename)
            arguments.append("%u %s %s" % (counter, process, filename))
            counter += 1
    
    print("Number of merging jobs: %u" % len(arguments))


    # Create jobdir and subdirectories
    jobdir = os.path.join(args["jobdir"], process,"merged")
    print("Jobdir: %s" % jobdir)
    mkdir(jobdir)
    mkdir(os.path.join(jobdir, "out"))
    mkdir(os.path.join(jobdir, "log"))
    mkdir(os.path.join(jobdir, "err"))

    # Write jdl file
    out = open(os.path.join(jobdir, "job_merged.jdl"), "w")
    out.write(jdl.format(PROCESS=process))
    out.close()

    # Write argument list
    arglist = open(os.path.join(jobdir, "arguments.txt"), "w")
    for a in arguments:
        arglist.write(a+"\n")
    arglist.close()

    # Write job file
    jobfile = open("job_merge.sh", "r").read()
    job = open(os.path.join(jobdir, "{PROCESS}.sh".format(PROCESS=process)), "w")
    job.write(jobfile)
    job.close()


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
