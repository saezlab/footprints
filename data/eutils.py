#!/usr/bin/env python2

# write scripts to check integrity of index.txt using die NCBI eutils and gds database

import os
import glob
import re
import pandas as pd
import pandas.io.parsers as parsers
import pickle
import ftplib
import subprocess
from Bio import Entrez
Entrez.email = "ms2126@cam.ac.uk"

# gets an NCBI record for a given GSE
def gse2gsm(gse, deep=False): #TODO: deep not implemented
    handle = Entrez.esearch(db="gds", term=gse+" AND GSM[Filter]", retmax=1000) # default: 20
    record = Entrez.read(handle)
    handle.close()
    idList = record['IdList']

    handle = Entrez.esummary(db="gds", id=",".join(idList))
    record = Entrez.read(handle)
    handle.close()
    return record

# converts an NCBI record to a dict structure GPL>GSM>fields
def record2index(record):
    index = {}
    for gpl in set([r['GPL'] for r in record]):
        index['GPL'+str(gpl)] = {r['Accession']:r for r in record if r['GPL']==gpl}
    return index

# get GSEs from index.txt, construct GSE/GSM index structure
def constructLookup(ifile='index.txt'):
    index = parsers.read_csv(ifile, sep="\t")
    return {gse:record2index(gse2gsm(gse)) for gse in set(index['GSE'])}

# check if index agrees with metadata from NCBI eutils
def checkIndex():
    pass

# download all available combinations between GSE and GPL
def downloadArrays(gses, datadir='../DATA/RAW/Transcriptomic/SPEED/CELS'):
    ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    for gse in gses:
        gsefile = gse + "_RAW.tar"
        destfile = os.path.join(datadir, gsefile)

        # download all GSEs not already downloaded
        if not os.path.isfile(destfile):
            try:
                url = "/".join(["geo/series", gse[0:-3]+"nnn", gse, 'suppl', gsefile])
                ftp.retrbinary('RETR %s' % url, open(destfile, 'wb').write)
            except:
                os.remove(destfile)

    ftp.quit()

#+download from HTTP, are there GSEs that are not on ftp?

def processArrays(lookup, method='rma', datadir='../DATA/RAW/Transcriptomic/SPEED/CELS'):
    rscript = os.path.join(os.getcwd(), 'CEL2RData.r')
    os.chdir(datadir)
    for gse,val in lookup.iteritems():
        fname = gse + "_RAW.tar"
        print("*** extracting " + fname)
        if os.path.isfile(fname) and not subprocess.call("tar -xf %s" % fname, shell=True):
            for gpl in val.keys():
                filenames = lookup[gse][gpl].keys()
                outfile = gse + "_" + gpl + "_" + method + ".RData"
                if not os.path.isfile(outfile):
                    subprocess.call("Rscript {} {} {} {} {}".format(
                        rscript, datadir, method, outfile, " ".join(filenames)), shell=True)

        subprocess.call("rm -f *.gz *.CEL", shell=True)

if __name__ == '__main__':
    lookup = constructLookup()
    downloadArrays(lookup.keys())
#    check()
    wd = os.getwd()
    processArrays(lookup, "rma") # TODO: fix dirchange
    os.chdir(wd)
    processArrays(lookup, "frma")
    os.chdir(wd)


# older code:
#######
def getGSMforGSE(gse):
    gsedict = {idList[i]:entry for i,entry in enumerate(record)}

    for myid in idList:
        handle = Entrez.efetch(db="gds", id=myid, retmode="text") # xml does not work here
        record = handle.read()
        handle.close()
        gsm = re.findall('(?<=Accession:\ )GSM[0-9]+', record)[0]
        gsedict[gsm] = gsedict.pop(myid)
        try:
            gsedict[gsm]['FTP'] = re.findall('ftp://[a-zA-Z0-9\/\.]+',record)[0]
        except:
            gsedict[gsm]['FTP'] = None

    return gsedict


if False: #__name__ == '__main__':
    # read index.txt in pandas dataframe (this should have headings already)
    index = parsers.read_csv('index.txt', sep="\t")
    allGSEs = set(index['GSE'])

    ##############
    # perform sanity checks listed below
    #############

    entries = {}
    try:
        with open('entries.pickle', 'rb') as fp:
            entries = pickle.load(fp)
        print("reading contents from pickled entries file")
    except:
        for gse in allGSEs:
            print("fetching metadata for {0}".format(gse))
            entries[gse] = getGSMforGSE(gse)
        with open('entries.pickle', 'wb') as fp:
            pickle.dump(entries, fp)

    # associate GSEs with GSMs
    gsegsm = {}
    for gse in allGSEs:
        gsegsm[gse] = entries[gse].keys()

    # check if any of the GSMs listed for a GSE do not actually belong to it
    for i,gse,gsm in index[['GSE','GSM']].itertuples():
        if gsm not in gsegsm[gse]:
            print("Warning: {0} not in {1}".format(gsm, gse))

    # check right platforms for all GSMs
    for i,gse,gsm,gpl in index[['GSE','GSM','GPL']].drop_duplicates().itertuples():
        gplentrez = 'GPL'+entries[gse][gsm]['GPL']
        if gpl != gplentrez:
            print("Warning: {0} {1}, {2} in .txt file but {3} in Entrez".format(gse,gsm,gpl,gplentrez))

    # check that in each subset controls are matched up with experiments
    for subset in set(index['SUBSET']):
        myset = index[index['SUBSET']==subset]

    #    eff = list(myset['EFFECT'])
    #    if eff.count('control') != eff.count('activating') \
    #            and eff.count('control') != eff.count('inhibiting') \
    #            and eff.count('control') != eff.count('unknown'): # tamoxifen treatments FIXME
    #        print("Warning: control/treatment not matched for subset {}".format(subset))

        same = myset[['GSE','GPL','SOFT_NUMSPACE','AVAIL','PATHWAY','CELLS','TREATMENT','TIME','COMMENT']]
        if len(same.drop_duplicates()) != 1:
            print("Warning: factors not constant for subset {}".format(subset))

    # platform always log or linear
    ## this is only relevant if processed data (eg. SOFT format)

    # check data availability for all GSMs
    ## implicitly did this with downloading; should maybe add it to index?

