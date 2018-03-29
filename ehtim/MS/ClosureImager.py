#!/usr/bin/env python
import numpy as np
import optparse
from DDFacet.Other import MyLogger
log = MyLogger.getLogger("ClosureImager")
import ClassWrapEHTImager
import pickle

def read_options():
    desc="""Run EHT imager"""
    
    opt = optparse.OptionParser(usage='Usage: %prog <options>',version='%prog version 1.0',description=desc)
    
    group = optparse.OptionGroup(opt, "* Data selection options")
    group.add_option('--MSName',type=str,help='',default=None)
    group.add_option('--ColName',type=str,help="",default='DATA')
    group.add_option('--FITSName',type=str,help="",default=None)
    group.add_option('--FlagAnts',type=str,help="Default is %default",default="CS,RS1,RS2,RS3,RS4")
    
    opt.add_option_group(group)
    options, arguments = opt.parse_args()
    f = open("last_param.obj","wb")
    pickle.dump(options,f)
    return options

        
#########################################

def main(options=None):
    if options is None:
        f = open("last_param.obj",'rb')
        options = pickle.load(f)
    
    MM=ClassWrapEHTImager.ClassWrapEHTImager(**options.__dict__)


if __name__=="__main__":
    OP=read_options()
    main(OP)



