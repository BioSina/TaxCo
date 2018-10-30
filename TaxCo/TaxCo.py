'''
TaxCo

@author: Sina Beier
'''

import argparse
import Convert
import Correlate
import Metadata
import Networking
import Plot
import Filter
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    ingroup = parser.add_argument_group(title='Input files (mutually exclusive)')
    group = ingroup.add_mutually_exclusive_group(required=True)
    group.add_argument("-qiime",metavar = ("<L2>","<L3>","<L4>","<L5>","<L6>"), type=str,nargs=5, help="Paths to QIIME taxon summaries, sorted from L2 to L6. Samples have to be in the same order in each file.")
    group.add_argument("-mothur",metavar="<tax.summary file>", type=str, nargs=1, help="Path to mothur tax.summary file (including all samples)")
    group.add_argument("-megan",metavar = ("<L2>","<L3>","<L4>","<L5>","<L6>"), type=str,nargs=5, help="Paths to MEGAN CSV files, sorted from L2 to L6. Collapse tree to matching level and export the summarized countsxs.")
    parser.add_argument("-out",metavar="<output directory>", type=str, default=os.getcwd()+"/TaxCoOut",help="Path to output directory")
    parser.add_argument("-pval",metavar="<p-value>", type=float,default=0.01, nargs=1, help="P-value cutoff")
    parser.add_argument("-cutoffs",metavar="<decimal number>", type=float, nargs='+' ,default=0.6,help="Cutoffs for correlation strength, with a dot as decimal separator")
    parser.add_argument("-coefficients", metavar="<correlation coefficients>", type=str, nargs='+', default="kendall", help="One or more correlation coefficients (kendall, spearman, pearson)")
    parser.add_argument("-metadata", metavar="<metadata>", type=str, help="Metadata in the format MEGAN exports it in (tab separated, samples in rows, data in columns)")
    
    args = parser.parse_args()
    
    #Convert
    if args.qiime:
        Convert.QIIMEConvert(args.qiime[0],args.qiime[1] ,args.qiime[2], args.qiime[3] ,args.qiime[4] , args.out)
    if args.mothur:
        Convert.mothurConvert(args.mothur[0], args.out)
    if args.megan:
        print("Conversion")
        Convert.MEGANConvert(args.megan[0],args.megan[1] ,args.megan[2], args.megan[3] ,args.megan[4] , args.out)
    #Correlate
    for i in args.coefficients:
        c = i.lower()
        print("Correlation on "+c)
        Correlate.correlate(args.out,c,args.pval,0.0)
        print("Metadata for "+c)
        Metadata.correlate(args.out,args.metadata,c,args.pval,0.0)
    print("Networking")
    Networking.writeNetworks(args.out)
    #Filter
    print("Filtering")
    Filter.filterAll(args.out,args.cutoffs)
    #Plot
    print("Plotting!")
    Plot.network(args.out)
    print("DONE!")