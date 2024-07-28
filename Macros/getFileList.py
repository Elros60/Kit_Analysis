#!/usr/bin/env python3

import argparse
import os
import alienpy.alien as alien
from glob import glob

"""
Macro for preparation of input text file used for AOD merger  
---- How to run:
python3 getFileList.py -r 505713 -y 2021 -p NOV -u arakotoz
"""

def getRunFileList(Input, Output, Local):

    if Local:
        #findCmd = "find {}/*/{}".format(Input, "AO2D.root")
        findCmd = "{}/**/{}".format(Input, "AO2D.root")
        print(findCmd)
        # FileList = os.listdir(findCmd)
        FileList = glob(findCmd, recursive=True)
        print(FileList)
        #results = out.ansdict['results']
        with open(Output, "a") as outFile:
            for dic in FileList:
                outFile.write("{}\n".format(dic))

    else:    
        alien.setup_logging()
        jal = alien.AliEn()
        findCmd = "find {} {} -e AOD".format(Input, "AO2D.root")
        print(findCmd)
        out = jal.run(findCmd)

        results = out.ansdict['results']
        with open(Output, "a") as outFile:
            for dic in results:
                outFile.write("alien://{}\n".format(dic["lfn"]))

        return out.exitcode



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get file lists")
    parser.add_argument("--Input", "-i", help="Directory to input files", required=True)
    parser.add_argument("--Output", "-o", help="Directory to output files", default="./input_data.txt")
    parser.add_argument("--Local", "-l", help="Local mode", default=True)

    args = parser.parse_args()
    getRunFileList(args.Input, args.Output, args.Local)
