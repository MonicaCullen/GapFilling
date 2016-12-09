#!/Usr/bin/python
# -*- coding:utf8 -*-
# @Date   :2016/10/14
# @Author :WuLing
# @Email  :wu_l@tib.cas.cn

import os 
import sys
import json
import copy
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import MolFromSmiles as mfsmi
from rdkit.Chem.AllChem import MolFromSmarts as mfsma
from rdkit.Chem.AllChem import MolToSmiles as mtsmi
from rdkit.Chem.AllChem import MolToSmarts as mtsma

sys.path.append('..')
from rxnpool import InitialiseNeutralisationReactions,NeutraliseCharges,DrawImage

def txt2json(filepath = False):
    
    if not filepath:
        filepath = './no_produce.txt'

    txtFile = open(filepath)

    keyValues = dict()
    keyList = list()
    valuesList = list()

    n = 0
    
    for line in txtFile:

        key = line.split('\t')[0].strip()
        value = line.split('\t')[1].strip()

        keyList.append(key)

        if value != '*' and bool(value) != False :
            try :
                value = NeutraliseCharges(mtsmi(mfsmi(value)))[0]
                keyValues[key] = value
            except:
                # print 'error value %s' % n ,value
                n += 1
                continue
            
            if value.count('.') != 0:
                values = value.split('.')
                for val in values:
                    valuesList.append(val)
            else:
                valuesList.append(value)

    keyList = list(set(keyList))
    valuesList = list(set(valuesList))
    keyValues = dedup(keyValues)

    savepath1 = filepath.rsplit('.',1)[0] + '-keyValues.json'
    savepath2 = filepath.rsplit('.',1)[0] + '-keyList.json'
    savepath3 = filepath.rsplit('.',1)[0] + '-valuesList.json'
    
    with open(savepath1,'w') as f1:
        json.dump(keyValues,f1,indent = 2)

    with open(savepath2,'w') as f2:
        json.dump(keyList,f2,indent = 2)

    with open(savepath3,'w') as f3:
        json.dump(valuesList,f3,indent = 2)

    return(keyValues,keyList,valuesList)

def allDeNosyn():
    
    mets_ex_NonProduce = txt2json('./mets_ex_NonProduce.txt')[0]
    noSynName = json.load(open('./noSyn.json'))
    noSyn = dict()
    
    mets_ex_NonProduce_noSyn = copy.deepcopy(mets_ex_NonProduce)
    
    for i in noSynName:
        
        if i in mets_ex_NonProduce.keys():
           
            mets_ex_NonProduce_noSyn.pop(i)
            noSyn[i] = mets_ex_NonProduce[i]


    with open('./mets_ex_NonProduce-noSyn.json','w') as fk:
        json.dump(mets_ex_NonProduce_noSyn,fk,indent = 2)
    with open('./noSyn-keyValues.json','w') as fn:
        json.dump(noSyn,fn,indent = 2)

    return (mets_ex_NonProduce_noSyn,noSyn)
   
def dedup(f):
    
    _f = dict()
    for i,j in f.items():
        if j not in _f.values():
            _f[i] = j
    print 'before deduplication',len(f)
    print 'After  deduplication',len(_f)
    return _f
def constant():

    C = txt2json('./mets_ex_NonProduce.txt')[2] #list
    A = txt2json('./no_produce.txt')[2] #list
    C_B,B = allDeNosyn() #dict

    allcp = A + C
    noPro = A + B.values()

    return(A,B,C,C_B,allcp,noPro)


def validPred(Choice = 1):
    '''
    mets_ex_NonProduce:C, NonProduce:A , noSyn:B, +:& , -:_
       mets_ex_NonProduce_noSyn :C_B
        A,C,B分别代表三个化合物的集合 '''

    choice = {1:'at least one ps in A+B',
              2:'at least one ps in Nopro and remained ps in C_B '}

    (A,B,C,C_B,allcp,noPro) = constant()
    
    resultDir = '../preResult/'
    f = open('./%s.txt' % choice[Choice],'w') 

    for filename in os.listdir(resultDir):

        for name in os.listdir(os.path.join(resultDir+filename)):

            path = resultDir+filename+'/'+name
            
            if name == 'result.txt':
                result_txt = open(path).read()
                query = mtsmi(mfsmi(result_txt.split('\n',1)[0].split(':',1)[1].strip()))
            
            if name == 'result.json':
                resultDict = json.load(open(path))

        _resultDict = copy.deepcopy(resultDict)
        
        for i,item in resultDict.items():
            smirks = item["Smikrs"]
            rs = smirks.split(">>")[0].split(".")
            ps = smirks.split(">>")[1].split(".")

            can = True

            for r in rs :
                if r == '[OH2]' or r == '[H]O[H]':
                    r = 'O'
               
                r = NeutraliseCharges(mtsmi(mfsmi(r)))[0]
                
                if r not in C_B.values():
                    can = False
                    _resultDict.pop(str(i))
                    break
                else:
                    pass

            if can:
                if Choice == 1:
                    can = False
                    for p in ps:
                        if p == '[OH2]' or p == '[H]O[H]' :
                            p = 'O'
                        try:
                            p = NeutraliseCharges(mtsmi(mfsmi(p)))[0]
                        except:
                            p = NeutraliseCharges(mtsmi(mfsma(p)))[0]
                        if p in noPro:
                            can = True
                        else:
                            pass
                    if not can:
                        _resultDict.pop(str(i))

                elif Choice == 2:
                    can = False
                    remained = True
                    for p in ps:
                        if p == '[OH2]' or p == '[H]O[H]' :
                            p = 'O'
                        try:
                            p = NeutraliseCharges(mtsmi(mfsmi(p)))[0]
                        except:
                            p = NeutraliseCharges(mtsmi(mfsma(p)))[0]
                        if p in noPro:
                            can = True
                        elif p in C_B.values():
                            pass
                        else:
                            remained = False
                    if (not can) or (not remained):
                        _resultDict.pop(str(i))

        output = filename+'\t'+'before:'+str(len(resultDict.keys()))+'\t'+'after:'+str(len(_resultDict))+'\n'
        
        f.write(output)
        f.flush()

        with open(os.path.join(resultDir,filename,'result_%s.json' % Choice),'w') as fn:
            json.dump(_resultDict,fn,indent = 2)

       # for j,_item in _result_json.items():
       #   rxnSmiFile = os.path.join(resultDir+filename+'/'+str(j)+".smi")
       #    rxnImagefile = os.path.join(resultDir+filename+'/'+str(j)+".png")
       #     DrawImage(rxnSmiFile,rxnImagefile) 
    f.close()

def validPredSum(Choice = 1):

    (A,B,C,C_B,allcp,noPro) = constant()

    resultDir = '../preResult/'

    preNoPro_querySmirks = dict()

    for filename in os.listdir(resultDir):

        for name in os.listdir(os.path.join(resultDir+filename)):

            path = resultDir+filename+'/'+name
            
            if name == 'result.txt':
                result_txt = open(path).read()
                query = mtsmi(mfsmi(result_txt.split('\n',1)[0].split(':',1)[1].strip()))
            
            if name == 'result_%s.json' % Choice: 
                resultDict = json.load(open(path))
    
        for i,item in resultDict.items():
            
            smirks = item["Smikrs"]
            rs = smirks.split(">>")[0].split(".")
            ps = smirks.split(">>")[1].split(".")

            for p in ps :

                if p == '[OH2]' or p == '[H]O[H]' :
                    p = 'O'
                try:
                    p = NeutraliseCharges(mtsmi(mfsmi(p)))[0]
                except:
                    p = NeutraliseCharges(mtsmi(mfsma(p)))[0]
                
                if p in noPro:
                    if p not in preNoPro_querySmirks.keys():
                        preNoPro_querySmirks[p] = dict()
                        preNoPro_querySmirks[p][query] = dict()
                        preNoPro_querySmirks[p][query][i] = smirks
                
                    else:
                        if query in preNoPro_querySmirks[p].keys():
                            preNoPro_querySmirks[p][query][i] = smirks
                        else:
                            preNoPro_querySmirks[p][query] = dict()
                            preNoPro_querySmirks[p][query][i] = smirks
                else:
                    pass

    with open('./preNoPro_querySmirks_%s.json' % Choice,'w') as fn:
        json.dump(preNoPro_querySmirks,fn,indent = 2)

     
def main():
    #efPredict3()
    IMAGEDIR = './preNoPro_querySmirks_2_PNG/'

    if not os.path.exists(IMAGEDIR):
        os.makedirs(IMAGEDIR)

    f = json.load(open("./preNoPro_querySmirks_1.json"))

    savefile = open("./preNoPro_querySmirks_1_summary.txt",'w')
   
    savefile.write('noPro ps types:' +'\t'+str(len(f.keys()))+'\n')
    savefile.write('-'*50 + '\n')
    m = 0
    num = 1
    for ps,item1 in f.items():
        savefile.write('ps%s: ' % num + ps +'\n')
        savefile.write('\t'+'predRs types: ' + str(len(item1.keys()))+'\n')
        num += 1

        n = 0
        for rs,item2 in item1.items():
            savefile.write('\t'*2 + 'smirks types: ' + str(len(item2.keys()))+'\n')
            n = n + len(item2.keys())

            # for _id,smirks in item2.items():
            #     with open(os.path.join(IMAGEDIR,'%s.smi' % _id),"w") as fk:
            #         fk.write(smirks.strip())
                
            #     rxnSmiFile = os.path.join(IMAGEDIR,'%s.smi' % _id)
            #     rxnImagefile = os.path.join(IMAGEDIR,"%s.png" % _id)
            #     DrawImage(rxnSmiFile,rxnImagefile)

        m += n
        
        savefile.write('total smirks: ' + str(n)+'\n')
        savefile.write('-'*50 + '\n')
    savefile.write('all ps total smirks: '+ str(m) + '\n')
if __name__ == '__main__':
    main()
    # validPred(Choice = 2)
    # validPredSum(Choice = 2)
