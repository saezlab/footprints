#!/usr/bin/env ipython2

def table2dict(table):
    mydict = { 'time':'undefined duration','comment':'none' }
    for line in table.strip().split('\n'):
        cc = line.strip().split('\t')
        mydict[cc[0].strip()] = cc[-1].strip()
    return mydict

import re
import os
index = 1

#directories = ! find . -type d -mindepth 2 -maxdepth 2

fields = ['ID', 'GSE', 'SUBSET', 'GSM', 'GPL', 'SOFT_NUMSPACE', 'AVAIL', 'FILE',
          'PATHWAY', 'CELLS', 'EFFECT', 'TREATMENT', 'TIME', 'COMMENT']

with open("index.txt", "w") as fp:
    fp.write("\t".join(fields))
    fp.write("\n")


for dd in directories:
    pathway = dd.split('/')[1]
    gse = re.findall("GSE[0-9]+", dd)[0]
    subset = re.findall("GSE[0-9a-zA-Z]+", dd)[0]

    with open(dd + "/Description.txt") as fp:
        desc = table2dict(fp.read())

    with open(dd + "/set.dat") as fp:
        res = fp.read().strip().split()
        if len(res) != 3:
            print(res,dd)
            raise AssertionError()
        gsetest,norm,gpl = [r.strip() for r in res]
        if gse != gsetest:
            print(pathway,gse,res,dd)
            raise AssertionError()

    with open(dd + "/gsm_accession.dat") as fp:
        lines = re.split('\r|\n', fp.read().strip())
        for line in lines:
            gsm,exp = line.strip().split()
            exptype = 'ERROR'
            if exp == '1':
                exptype = 'control'
            elif exp == '2':
                if 'effect' not in desc.keys():
                    print dd
                    raise AssertionError()
                exptype = desc['effect']
            else:
                raise AssertionError()

            log = "log"
            normalisation = "raw"
            filename = os.path.join("CELS", gse, gsm+".CEL")
            if not os.path.isfile(os.path.join("..", filename)):
                filename = "NA"
            if norm == "log": # is this true?
                normalisation = "raw"
                log = "linear" # in speed, log is a flag that we should log the data, ie. it is in linear space
            if norm == "normalized":
                normalisation = "normalized"

            with open("index.txt", "a") as fp:
                mylist = ["ARR%05i"%index,gse,subset,gsm,gpl,log,normalisation,filename,pathway,desc['cell line'],exptype,desc['modification'],desc['time'],desc['comment']]
                fp.write("\t".join([x.strip() for x in mylist]))
                fp.write("\n")
                index+=1

