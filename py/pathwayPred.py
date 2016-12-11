# !/usr/bin/python
#  _*_ coding:utf-8 _*_
# 2016/12/06
# Author:LingWu
# Email:wu_l@tib.cas.cn 

import os
import json
from crawler import *

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

def bridgePre(n):

    cpAddRxn = json.load(open('../modules/iJO1366/pathwayVerify/cpAddRxn.json'))
    keggIdModelId = json.load(open('../modules/iJO1366/original_model/KeggIdModelId.json'))
    metsAllSmi = json.load(open('../modules/iJO1366/original_model/iJO1366_metsAll_smi(kegg18).json'))

    cpAddRxn_Smi = dict()

    for cp,addRxns in cpAddRxn.items():
        print cp

        if len(addRxns) == n  :
            try:
                cpModelId = keggIdModelId[cp]
                cpSmi =  metsAllSmi[cpModelId]
            except:
                space = 1
                moltxt = crawlerKEGGCP(cp)

                head = moltxt.split(' ')[0].strip()
                if len(head) == 1:
                    space = 2

                f1 = open('./tmp.mol' ,'w')
                f1.write('\n\n\n'+' '*space)
                f1.write(moltxt)
                f1.flush()
                f1.close()
                mol = AllChem.MolFromMolFile('./tmp.mol')
                cpSmi = AllChem.MolToSmiles(mol)
            
            if cpSmi:
                cpAddRxn_Smi[cpSmi] = list()
                
                for rxn in addRxns:

                    rxnSmi = dict()
                    for p,coef in rxn.items():
                        
                        if p.startswith('add'):
                            pKeggId = p.split('add_met_kegg_')[1].split('_')[0]
                            if pKeggId == cp:
                                pSmi = cpSmi
                                rxnSmi[pSmi] = coef
                                continue
                            else:
                                try:
                                    pModelId = keggIdModelId[pKeggId]
                                    pSmi =  metsAllSmi[pModelId]
                                except:
                                    space = 1
                                    moltxt = crawlerKEGGCP(pKeggId)

                                    head = moltxt.split(' ')[0].strip()
                                    if len(head) == 1:
                                        space = 2

                                    f1 = open('./tmp.mol' ,'w')
                                    f1.write('\n\n\n'+' '*space)
                                    f1.write(moltxt)
                                    f1.flush()
                                    f1.close()
                                    mol = AllChem.MolFromMolFile('./tmp.mol')
                                    pSmi = AllChem.MolToSmiles(mol)
                                    rxnSmi[pSmi] = coef

                        else:
                            pSmi = metsAllSmi[p]
                            rxnSmi[pSmi] = coef
                    if len(rxnSmi.keys()) == len(rxn.keys()):
                        cpAddRxn_Smi[cpSmi].append(rxnSmi)
                    else:
                        print ' %s rxn %s mets in rxn have no smi' % (cp,str(addRxns))

    with open('../modules/iJO1366/pathwayVerify/cpAddRxn_Smi.json', 'w') as fn:
        json.dump(cpAddRxn_Smi,fn,indent = 2)
    print len(cpAddRxn_Smi.keys())
    
def main():
    # getExtendRxn()
    # bridgePre(1)

    # f = json.load(open('../modules/iJO1366/pathwayVerify/cpAddRxn_Smi.json'))
    # print len(f.keys())
if __name__ == '__main__':
    main()
