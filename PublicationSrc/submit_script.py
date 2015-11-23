#!/usr/bin/python                                                                                                          
import os
import sys

JOBDIR="/ebio/ag-neher/home/rneher/Projects/Dan_Balick/Ratchet/PublicationSrc/"
command = "python "+ JOBDIR+"ratchet.py "
for arg in sys.argv[1:]:
    command+=" "+arg

print command
os.system(command)


