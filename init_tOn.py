#!/usr/bin/python2

import shlex
import numpy as np
import sys
import subprocess 

def get_exitcode_stdout_stderr(cmd):
    """
    Execute the external command and get its exitcode, stdout and stderr.
    """
    args = shlex.split(cmd)

    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    exitcode = proc.returncode
    #
    return exitcode, out, err

path = sys.argv[1]
spath= sys.argv[2]
_skipper=0
_initvals=	(np.genfromtxt(path+"/x", skip_header =_skipper, delimiter="\t", dtype='float')[:]) #_num.. makes the var internal == private
_initparams=	(np.genfromtxt(path+"/parameter", skip_header =_skipper, delimiter=None, dtype='float')[:]) #_num.. makes the var internal == private
jstart=_initparams[21]*1E6
jstop=jstart+10
feed_P=_initparams[22]
feed_A=_initparams[23]
histOm=_initparams[-1]

t_stop=100

#
argsstring="-spath "+ spath+" -Einit_Re "+str(_initvals[0])+" -Einit_Im " + str(_initvals[1]) + " -Rhoinit " + str(_initvals[2]) + " -Ninit " + str(_initvals[3]) +" -rmlog"+" -jstart " + str(jstart) + " -jstop "+str(jstop) + " -fA_ss " + str(feed_A) + " -fP_ss " +str(feed_P)+ " -historyomega " + str(histOm) + " -t_stop " + str(t_stop) + " -nolog" 
#argsstring="-spath "+ spath+" -rmlog"+" -jstart " + str(jstart) + " -jstop "+str(jstop)+ " -dt " +str(0.001)
cmd="./testben "+argsstring
print cmd

#subprocess.call(cmd, shell=True)

exitcode, out, err= get_exitcode_stdout_stderr(cmd)
print out, err


