# !/usr/bin/python
#  _*_ coding:utf-8 _*_
# 2016/12/06
# Author:LingWu
# Email:wu_l@tib.cas.cn 

import os
import json
def getExtendRxn():

    fileDir = '../modules/iJO1366/pathwayVerify/model_id_changed/'

    cpAddRxn = dict()
    for filename in os.listdir(fileDir):
        cpId = filename.split('_')[0]
        filepath = os.path.join(fileDir,filename)
        rxns = json.load(open(filepath))

        for rxn in rxns["reactions"]:
            if rxn["id"].count('add_rxn'):
                if cpId not in cpAddRxn.keys():
                    cpAddRxn[cpId] = list()

                metabolites = rxn['metabolites']
                cpAddRxn[cpId].append(metabolites)
    with open(os.path.join('../modules/iJO1366/pathwayVerify/','cpAddRxn.json'),'w') as fn:
        json.dump(cpAddRxn,fn,indent = 2)
    print cpAddRxn

def main():
    # getExtendRxn()

    cpAddRxn = json.load(open('../modules/iJO1366/pathwayVerify/cpAddRxn.json'))
    n = 0
    print len(cpAddRxn.keys())
    for cp,addRxns in cpAddRxn.items():
        if len(addRxns) ==14 :
            n += 1
    print n

    '''
    total:257; 
    n = 1,77; 
    n = 2,53; 
    n = 3,38;  
    n = 4,24, 
    n = 5,13; 
    n = 6,11; 
    n = 7,10;
    n = 8,11; 
    n = 9,5; 
    n = 10,3; 
    n = 11,5;
    n = 12,4;
    n = 14,1;
    n = 15,1;
    n = 16,1;
    '''
if __name__ == '__main__':
    main()
