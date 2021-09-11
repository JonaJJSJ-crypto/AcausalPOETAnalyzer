#!/usr/bin/env python

import os
import sys

THEEOSDIR = "/eos/home-e/ecarrera/lee-wick/signal/"
NUMJOBS = 1000
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
    for i in range(0,NUMJOBS):
        genfile = GENTYPE+"_"+theprocess+"_"+str(i)+".root"
        hltfile = HLTTYPE+"_"+theprocess+"_"+str(i)+".root"
        recofile = RECOTYPE+"_"+theprocess+"_"+str(i)+".root"
        if os.path.exists(thedir+genfile):
            print(thedir+genfile+"\t exists." )
        else:
            print(thedir+genfile+"\t DOES NOT EXIST!!!!!!" )
        if os.path.exists(thedir+hltfile):
            print(thedir+hltfile+"\t exists." )
        else:
            print(thedir+hltfile+"\t DOES NOT EXIST!!!!!!" )
        if os.path.exists(thedir+recofile):
            print(thedir+recofile+"\t exists." )
        else:
            print(thedir+recofile+"\t DOES NOT EXIST!!!!!!" )

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
