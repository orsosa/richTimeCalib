#!/usr/bin/env python
import glob
import os
from sys import argv
import subprocess

RICHTIMECALIB = os.getenv("RICHTIMECALIB")
if RICHTIMECALIB == "":
    print "environment variable $RICHTIMECALIB must be set to"
    print "the location of the time calibration suite."
    exit(1)

############## global setting ###########
outdir_pref = "/volatile/clas/claseg2/osoto/RICH/timeCalib/Hists"
appHist   = RICHTIMECALIB + "/bin/richTiming_app"
aBin = [RICHTIMECALIB + "/bin/richTimeOffsets"]
aBin.append(RICHTIMECALIB + "/bin/richTimeWalks")
aBin.append(RICHTIMECALIB + "/bin/richTimeCheck")
scriptHist= RICHTIMECALIB + "/script/run_richTiming.sh"
aScript  = [RICHTIMECALIB + "/script/run_TO_corr.sh"]
aScript.append(RICHTIMECALIB + "/script/run_TW_corr.sh")
aScript.append(RICHTIMECALIB + "/script/run_TC.sh")
maps = [RICHTIMECALIB + "/maps/SspRich_mapCHANNEL2PIXEL.txt"]
maps.append(RICHTIMECALIB + "/maps/SspRich_mapFIBER2PMT_sortbyPMT.tx")
ofile_pref = "RichTimeCalib_"
WF    = "rich_timeCalib"
DEBUG = False
RN    = 0
MAXJOBS = 2
TRACK = "debug"
#################################################
def checkdir(path):
    if not os.path.isdir(path):
        print "The directory " + path + " does not exist"
        print "creating it"
        cmd = "mkdir -p " + path
        subprocess.call(cmd,shell=True)
        return 0
        
####################

def get_file_list(path,regexp="*.hipo"):
    flist=glob.glob(path + "/" + regexp)
    return flist

###### import file list from file ############
def import_flist(fname):
    infile=open(fname)
    flist=[]
    for line in infile:
        flist.append(line.strip())
    return flist

##### add histogram job to workflow
def add_hist_job(wf,fname,phase=0,c=0):
    global outdir_pref, appHist, scriptHist, RN
    script = scriptHist
    jname = wf + "_" + "ph" + str(phase) + "_" + str(c)
    size = os.path.getsize(fname)/1024./1024.*1.15
    cmd  = "swif add-job -workflow " + wf + " -ram 3000mb -project clas12 -time 10min -disk " + "{0:.0f}".format(size) + "mb "
    cmd += " -track "+ TRACK
    cmd += " -phase " + str(phase)
    cmd += " -name " + jname
    cmd += " -shell /bin/bash" 
    cmd += " -input " + appHist.split("/")[-1] + " file:" + appHist
    cmd += " -input " + script.split("/")[-1] + " file:" + script
    cmd += " -input " + fname.split("/")[-1] + " file:" + fname
    T = int(phase/2)
    corr_type = ["Offsets","Walks"]
    if phase > 0:
        cmd += " -input richTime" + corr_type[T-1] + "_" + RN + ".out" + " file:" + outdir_pref + "/T" + corr_type[T-1] + "/richTime" + corr_type[T-1] + "_" + RN + ".out"

    outdir = outdir_pref + "/T" + str(T)
    checkdir(outdir)

    cmd += " ./" + script.split("/")[-1] + " " + str(T) +  " " + outdir + " " + RN
    if DEBUG : print (cmd)
    subprocess.call(cmd,shell=True)

#### add analysis job ##########################
def add_ana_job(wf,flist,phase=0):
    global outdir_pref, aBin, aScript, RN
    size = 24*4.
    cmd  = "swif add-job -workflow " + wf + " -ram 500mb -project clas12 -time 20min -disk " + "{0:.0f}".format(size) + "mb "
    cmd += " -track "+ TRACK
    cmd += " -phase " + str(phase)
    cmd += " -shell /bin/bash" 

    T = int((phase-1)/2)
    script = aScript[T]
    app = aBin[T]

    cmd += " -input " + app.split("/")[-1] + " file:" + app
    cmd += " -input " + script.split("/")[-1] + " file:" + script

    jname = wf + "_" + "ph" + str(phase) 

    corr_type = ["TOffsets","TWalks","TCheck"]
    outdir = outdir_pref + "/" + corr_type[T]
    jname += "_" + corr_type[T]

    cmd += " -name " + jname

    checkdir(outdir)

    c = 0
    for fname in flist:
        fname = fname.split("/")[-1]
        fname = outdir_pref + "/T"+ str(T) +"/" + fname.replace(".hipo","") + "__RichTimeCalibE_" + RN + "_" + str(T) +".root"
        cmd += " -input " + fname.split("/")[-1] + " file:" + fname
        c +=1
        if c>=MAXJOBS: break

    cmd += ' ./' + script.split("/")[-1] + " " + outdir + " " + RN
    if DEBUG : print (cmd)
    subprocess.call(cmd,shell=True)

######### create workflow #####
def create_wf():
    global WF
    cmd="swif create " + WF
    if DEBUG : print (cmd)
    subprocess.call(cmd,shell=True)

#### run workflow ############
def run_wf():
    global WF
    cmd="swif run -workflow " + WF
    subprocess.check_call(cmd.split(" "))

######## main routine ####

def main():
    global WF, RN, outdir_pref, MAXJOBS
    if len(argv)<2: 
        print "A directory with calibration files or filelist must be supplied"
        exit (1)
    if not os.path.exists(argv[1]):
        print "The file or path " + argv[1] + "does not exist"
        exit (1)

    if len(argv) == 3:
        try:
            RN = int(argv[2])
        except ValueError:
            print "The run number must be an integer!. You supplied " + argv[2]
            exit(0)

    if RN == 0:
        RN = "" 
    else:
        RN = str(RN)


    outdir_pref += "/" + RN
    create_wf()
    flist = []
    if os.path.isdir(argv[1]):
        flist = get_file_list(argv[1])
    elif os.path.isfile(argv[1]):
        flist = import_flist(argv[1])
    else:
        print argv[1] + " is not a file nor a directory."

    ##### Setting jobs phase ###########
    phase = 0
    ##### Putting jobs to make histograms without correction ####
    c=0
    for fn in flist:
        add_hist_job(WF,fn,phase,c)
        c +=1
        if c>=MAXJOBS: break

    ##### Setting jobs phase ###########
    phase +=1
    #### Putting jobs to make TO corrections ####
    add_ana_job(WF,flist,phase)

    ##### Setting jobs phase ###########
    phase +=1
    #### Putting jobs to make histogram with TO corrections ####
    c=0
    for fn in flist:
        print fn
        add_hist_job(WF,fn,phase,c)
        c +=1
        if c>=MAXJOBS: break

    ##### Setting jobs phase ###########
    phase +=1
    #### Putting jobs to make TW corrections ####
    add_ana_job(WF,flist,phase)

    ##### Setting jobs phase ###########
    phase +=1
    #### Putting jobs to make histogram with TO and TW corrections ####
    c=0
    for fn in flist:
        print fn
        add_hist_job(WF,fn,phase,c)
        c +=1
        if c>=MAXJOBS: break

    ##### Setting jobs phase ###########
    phase +=1
    #### Putting jobs to produce checking plots ####
    add_ana_job(WF,flist,phase)

if __name__=="__main__":
    main()
    run_wf()