#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
2018-03-08: modified default minPts to 3
2018-03-13: mode added for pre-set parameters
2018-03-26: modified cut option , removed
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"
__version__ = "0.9"

#sys library
import os, time, sys, logging, gzip, argparse

#glob settings
#epilog for argparse
EPILOG = "Any bug is welcome reported to caoyaqiang@picb.ac.cn, aidaosheng@picb.ac.cn, chenzhaoxiong@picb.ac.cn"


def getLogger(fn=None):
    """
    Creat the logger system.
    """
    #get the current time
    date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
    #if fn == None:
    #    fn = os.getcwd() + "/" + date.strip() + ".log"
    #set up logging, both write log info to console and log file
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(name)-6s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=fn,
        filemode='a')
    logger = logging.getLogger()
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.NOTSET)
    return logger


def callSys(cmds, logger=None):
    """
    Call systematic commands without return.
    """
    for c in cmds:
        try:
            logger.info(c)
        except:
            print c
        try:
            os.system(c)
        except:
            try:
                logger.error(c)
            except:
                print "ERROR!", c


def cFlush(r):
    """
    One line flush.
    """
    sys.stdout.write("\r%s" % r)
    sys.stdout.flush()


def mainHelp():
    """
    Create the command line interface of the main programme for calling loops.
    """
    epilog = EPILOG
    description = """
        Intra-chromosomal loops calling for ChIA-PET,HiChIP and high-resolution Hi-C data.
        """
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-f",
        dest="fnIn",
        required=True,
        type=str,
        help=
        "Input file is mapped PETs, BEDPE format in gzip. Replicates could be input as -f A.bedpe.gz,B.bedpe.gz,C.bedpe.gz. Loops will be called with pooled data."
    )
    parser.add_argument(
        "-o", dest="fnOut", required=True, type=str, help="Output prefix.")
    parser.add_argument(
        "-m",
        dest="mode",
        required=False,
        type=int,
        default=0,
        choices=[0, 1, 2, 3, 4],
        help=
        "Pre-set parameters and signicicance cutoff for different types of data. Default is 0, using values from -eps and -minPts. Set 1 for sharp peak like ChIA-PET data (CTCF,RAD21,eg..),set 2 for broad peak like ChIA-PET data (H3K27ac,H3K4me1 eg..), and set 3 for deep sequenced Hi-C (>200 million cis PETs), set 4 for HiChIP (>100 million cis PETs). Detail parameters will be logged in the log file."
    )
    parser.add_argument(
        "-eps",
        dest="eps",
        default=0,
        required=False,
        help=
        "Distance that define two points being neighbors, eps in cDBSCAN as key parameter. For sharp peak like ChIA-PET data (CTCF), it can be set as 1000,2000. For broad peak like ChIA-PET data, such as H3K27ac/H3K4me1, set it to 2000,5000. For data like HiChIP and Hi-C, set it larger,set several eps like 5000,7500,10000. Default is 0, cLoops can auto estimate a eps for initial result, maybe not good."
    )
    parser.add_argument(
        "-minPts",
        dest="minPts",
        default=0,
        help=
        "Points required in a cluster, minPts in cDBSCAN, key parameter. Empirically 5 is good for TFs and histone modification ChIA-PET data. For data like HiChIP and Hi-C, set it larger, like >=20. Since v0.9, it can be a seires, and the final loops will have the PETs>= max(minPts). For Hi-C data with ~200 million intra-chromosomal PETs, we set it to 20,30,40,50. For cohesin HiChIP data with ~30-40 million intra-chromosomal PETs, we set it to 10,15,20. You can custome it."
    )
    parser.add_argument(
        "-p",
        dest="cpu",
        required=False,
        default=1,
        type=int,
        help=
        "CPU number used to run the job, default is 1,set -1 to use all cpus available. Too many CPU could cause memory error."
    )
    parser.add_argument(
        "-c",
        dest="chroms",
        required=False,
        default="",
        type=str,
        help=
        "Whether to process limited chroms, specify it as chr1,chr2,chr3, default is not. Use this to filter reads in like chr22_KI270876v1_alt"
    )
    parser.add_argument(
        "-w",
        dest="washU",
        required=False,
        action="store_true",
        help=
        "Whether to save tracks of loops to visualize in washU. Default is No, set this flag to save."
    )
    parser.add_argument(
        "-j",
        dest="juice",
        required=False,
        action="store_true",
        help=
        "Whether to convert loops to 2d feature annotations to visualize in Juicebox. Default is No, set this flag to save."
    )
    parser.add_argument(
        "-s",
        dest="tmp",
        required=False,
        action="store_true",
        help=
        "Whether or not to save temp files for each chromosomes during processing. Set this flag for following calling differentially enriched loops or converting PETs to washU track or hic file load into juicebox. Default is not."
    )
    parser.add_argument(
        "-hic",
        dest="hic",
        required=False,
        action="store_true",
        help=
        "If input is HiChIP or high resolution Hi-C data, set this flag, using different significance cutoffs for loops than ChIA-PET data."
    )
    parser.add_argument(
        "-cut",
        dest="cut",
        required=False,
        default=0,
        type=int,
        help=
        "Initial distance cutoff to filter PETs, default is 0, only used for debuging."
    )
    parser.add_argument(
        "-plot",
        dest="plot",
        required=False,
        action="store_true",
        help=
        "Whether to plot estimated inter-ligation and self-ligation PETs distance distrbution, default is not."
    )
 
    parser.add_argument(
        "-v",
        dest="version",
        action="version",
        version="cLoops v%s" % __version__, )

    op = parser.parse_args()
    return op


def deloopHelp():
    """
    Create the command line interface for the script of deLoops 
    """
    epilog = EPILOG
    description = """
        Differentially enriched loops calling based on loops called by cLoops.
        For example:
        deLoops -fa a.loop -fb b.loop -da A -db B -p 10 
        """
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-fa",
        dest="fa",
        required=True,
        type=str,
        help=
        "Loops file called by cLoops. Only using significant loops as mark 1, you can change this in the .loop file."
    )
    parser.add_argument(
        "-fb",
        dest="fb",
        required=True,
        type=str,
        help="Loops file called by cLoops.")
    parser.add_argument(
        "-da",
        dest="da",
        required=True,
        type=str,
        help=
        "Directory for .jd files of loop file a, generated by cLoops with option -s 1."
    )
    parser.add_argument(
        "-db",
        dest="db",
        required=True,
        type=str,
        help=
        "Directory for .jd files of loop file b, generated by cLoops with option -s 1."
    )
    parser.add_argument(
        "-p",
        dest="cpu",
        required=False,
        default=1,
        type=int,
        help=
        "CPU number used to run the job, default is 1,set -1 to use all cpus available. Too many CPU could cause memory error."
    )
    parser.add_argument(
        "-c",
        dest="chroms",
        required=False,
        default="",
        type=str,
        help=
        "Whether to process limited chroms, specify it as chr1,chr2,chr3, default is not. Set it to the same one for cLoops."
    )
    parser.add_argument(
        "-dis",
        dest="dis",
        required=False,
        default=0,
        type=int,
        help=
        "Set a distance cutoff to filter PETs, could be the inter-ligation and self-ligation cutoff, default is 0."
    )

    op = parser.parse_args()
    return op


def jd2washUHelp():
    """
    Create the command line interface for the script of jd2washU.
    """
    epilog = EPILOG
    description = """
        Convert PETs level data to washU browser track for visualization. bedtools,bgzip,tabix are required.
        Example:
        jd2washU -d CTCF_ChIA-PET -o CTCF_ChIA-PET 
        """
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-d",
        dest="dir",
        required=True,
        type=str,
        help="Directory for .jd files, generated by cLoops with option -s 1.")
    parser.add_argument(
        "-o", dest="output", required=True, type=str, help="Output prefix.")
    parser.add_argument(
        "-ext",
        dest="ext",
        type=int,
        default=75,
        help=
        "Extension from the middle center of the PET to both ends,default is 75."
    )
    parser.add_argument(
        "-cut",
        dest="cut",
        type=int,
        default=0,
        help="Distance cutoff for PETs to filter, default is 0.")

    op = parser.parse_args()
    return op


def jd2juiceHelp():
    """
    Create the command line interface for the script of jd2juice.
    """
    epilog = EPILOG
    description = """
        Convert PETs level data to .hic file to load in juicebox. The command "juicer_tools pre" is required in the enviroment.
        For example:
        jd2juice -d CTCF -o test -org hg38
        """
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-d",
        dest="dir",
        required=True,
        type=str,
        help="Directory for .jd files, generated by cLoops with option -s 1.")
    parser.add_argument(
        "-o", dest="output", required=True, type=str, help="Output prefix.")
    parser.add_argument(
        "-org",
        dest="org",
        required=True,
        type=str,
        default="hg38",
        help="Organism required to generate .hic file,default is hg38.")
    parser.add_argument(
        "-cut",
        dest="cut",
        type=int,
        default=0,
        help="Distance cutoff for PETs to filter, default is 0.")
    op = parser.parse_args()
    return op


def jd2saturationHelp():
    """
    Create the command line interface for the script of jd2saturation.
    """
    epilog = EPILOG
    description = """
        Re-sampling PETs to estimate loops detection saturation. 
        For example:
        jd2saturation -jd CTCF_chr21/chr21-chr21.jd -o test -s 5 -eps 750 -minPts 5
        """
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-jd",
        dest="jd",
        required=True,
        type=str,
        help=
        "The .jd file used to re-sampling, generated by cLoops with option -s 1. Small chromosome is recommended,like chr21 for human or chr19 for mouse."
    )
    parser.add_argument(
        "-o", dest="output", required=True, type=str, help="Output prefix.")
    parser.add_argument(
        "-r",
        dest="repeats",
        type=int,
        default=5,
        help="Re-sampling times for each sequencing depth,default is 5.")
    parser.add_argument(
        "-s",
        dest="step",
        type=int,
        default=10,
        help="Re-sampling step,default is 10. ")
    parser.add_argument(
        "-eps",
        dest="eps",
        required=True,
        help="Eps to run. Eps will be the same for different sequencing depth. "
    )
    parser.add_argument(
        "-minPts",
        dest="minPts",
        required=True,
        type=int,
        help=
        "minPts to run. minPts will be re-caculated according to re-sampling sequencing depth, say 0.5 PETs and it equals 0.5*minPts."
    )
    parser.add_argument(
        "-cut",
        dest="cut",
        type=int,
        default=0,
        help="Distance cutoff for PETs to filter, default is 0.")
    parser.add_argument(
        "-p",
        dest="cpu",
        required=False,
        default=1,
        type=int,
        help=
        "CPU number used to run the job, default is 1,set -1 to use all cpus available. Too many CPU could cause memory error."
    )
    parser.add_argument(
        "-hic",
        dest="hic",
        required=False,
        default=False,
        type=bool,
        help=
        "If input is HiChIP or high resolution Hi-C data, set this to 1, using different significance cutoffs for loops than ChIA-PET data."
    )
    op = parser.parse_args()
    return op


def jd2fingerprintHelp():
    """
    Create the command line interface for the script of jd2saturation.
    """
    epilog = EPILOG
    description = """
        Get the finger print for the datasets using contact matrix with specific bin size, small bin sizes like 1000, 2000 are recommended.
        For exmaple:
        jd2fingerprint -jd CTCF_ChIA-PET,cohesin_HiChIP,HiC -o test -bs 2000 -plot 1 -p 10 -labels CTCF_ChIA-PET,cohesin_HiChIP,HiC
        """
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-d",
        dest="d",
        required=True,
        type=str,
        help=
        "The directory of cis .jd file used to draw quality control fingerprint,created by cLoops with option -s 1. Mutiple samples/datasets for comparasion should be assigned as -d a,b,c"
    )
    parser.add_argument(
        "-o", dest="output", required=True, type=str, help="Output prefix.")
    parser.add_argument(
        "-bs",
        dest="binSize",
        default=2000,
        type=int,
        help="Bin sizes for contact matrix, default is 2000.")
    parser.add_argument(
        "-labels",
        dest="labels",
        default="",
        help=
        "Labels for the datasets for ploting, default is the directory name.")
    parser.add_argument(
        "-plot",
        dest="plot",
        default=0,
        type=bool,
        help=
        "Set 1 for ploting the finger print, set 0 for not, default is not. ")
    parser.add_argument(
        "-p",
        dest="cpu",
        required=False,
        default=1,
        type=int,
        help=
        "CPU number used to run the job, default is 1,set -1 to use all cpus available."
    )
    op = parser.parse_args()
    return op
