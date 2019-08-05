#!/bin/bash
username=$(whoami)
USERDIR=/afs/cern.ch/user/${username:0:1}/$username
WORKDIR=/afs/cern.ch/work/${username:0:1}/$username

source ${USERDIR}/.bash_profile
source ${USERDIR}/AMS/install/amsvar_icc.sh

${WORKDIR}/Anomalous/anomalous $1 $2 # run!
