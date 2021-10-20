#!/usr/bin/env python

import os
import sys

THEEOSDIR = "/eos/home-e/ecarrera/lee-wick/signal/"
THEJOBDIR = "joblaunch"
NUMJOBS = 15000
GENTYPE = "gensimLW"
HLTTYPE = "hltLW"
RECOTYPE = "recoLW"


def parse_arguments():
    if not len(sys.argv) == 2:
        raise Exception("Run with './checkFiles.py <mass>'")
    return sys.argv[1]


def main(themass):
    theprocess = "LWSM"+themass+"DnR"
    thedir = THEEOSDIR+theprocess+"/" 
    print ("Testing directory: "+thedir)
    failedjobs = []
    for i in range(0,NUMJOBS):
        genfile = GENTYPE+"_"+theprocess+"_"+str(i)+".root"
        hltfile = HLTTYPE+"_"+theprocess+"_"+str(i)+".root"
        recofile = RECOTYPE+"_"+theprocess+"_"+str(i)+".root"
#        if os.path.exists(thedir+genfile):
#            print(thedir+genfile+"\t exists." )
#        else:
#            print(thedir+genfile+"\t DOES NOT EXIST!!!!!!" )
#        if os.path.exists(thedir+hltfile):
#            print(thedir+hltfile+"\t exists." )
#        else:
#            print(thedir+hltfile+"\t DOES NOT EXIST!!!!!!" )
        if not os.path.exists(thedir+recofile):
#            print(thedir+recofile+"\t DOES NOT EXIST!!!!!!" )
            failedjobs.append(i)

    #parse the arguments file in the joblaunch dir and create a new one with failed jobs for
    #manual relaunch
    print ("A total of "+str(len(failedjobs))+" have failed.  Do not worry, we will generate a failedarguments.txt file for resubmission.")
    thejobdir = THEJOBDIR+"/"+theprocess+"/"
    argsfile = "arguments.txt"
    theargsfile = thejobdir+argsfile
    failedargsfile = "failedarguments.txt"
    thefailedargsfile = thejobdir+failedargsfile
    ffargs = open(thefailedargsfile,"w")
    fargs = open(theargsfile,"r")
    for line in fargs.readlines():
        if (int(line.split()[0]) in failedjobs):
            ffargs.write(line)

    ffargs.close()

    #also create a new failedjob.jdl
    thejdlfile =  thejobdir+"job.jdl"
    thefailedjdlfile = thejobdir+"failedjob.jdl"
    os.system("cp "+thejdlfile+" "+thefailedjdlfile)
    thesed = 'sed -i -e "s,'+argsfile+','+failedargsfile+',g" '+thefailedjdlfile
    os.system(thesed)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
