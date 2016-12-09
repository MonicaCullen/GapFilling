# !/usr/bin/python
#  _*_ coding:utf-8 _*_
# 2016/11/01
# Author:LingWu
# Email:wu_l@tib.cas.cn

import os
import copy
import json
import cobra
import rdkit
from testV5 import *
from rdkit import Chem
from rdkit.Chem import AllChem
from cobra import Model, Reaction, Metabolite
from rdkit.Chem.AllChem import MolFromSmiles as mfsmi
from rdkit.Chem.AllChem import MolFromSmarts as mfsma
from rdkit.Chem.AllChem import MolToSmiles as mtsmi
from rdkit.Chem.AllChem import MolToSmarts as mtsma

def genUniversal():
    
    keggRxnfile_ID = open('./keggDatebase/rxnCompID_set.keg')
    keggRxnfile_NAME_1 = open('./keggDatebase/rxnCompName.keg')
    keggRxnfile_NAME_2= open('./keggDatebase/RxnIdEqName.txt')
    KeggIdModelId = getKeggIdModelId()

    rxnIdComId,allComId = tmp(keggRxnfile_ID,coefficient = True)
    rxnIdComName_1,allComName_1 = tmp(keggRxnfile_NAME_1,coefficient = False)
    rxnIdComName_2,allComName_2 = tmp(keggRxnfile_NAME_2,coefficient = False)

    keggRxnfile_ID.close()
    keggRxnfile_ID = open('./keggDatebase/rxnCompID_set.keg')


    print len(rxnIdComId),len(allComId)
    print len(rxnIdComName_1),len(allComName_1)
    print len(rxnIdComName_2),len(allComName_2)


    idName = dict()
    for rxnid,ComId in rxnIdComId.items():
        for i,com in enumerate(ComId):
            try:
                name = rxnIdComName_2[rxnid][i]
            except:
                name = rxnIdComName_1[rxnid][i]

            if name.count(' ') and name[:name.find(' ',1)].isalnum():
                name = name.split(' ',1)[1]
            idName[com] = name
    
    universal = Model('Universal_Reactions')
    
    reactions = list()
    
    for line in keggRxnfile_ID:

        rxnId = line.split('\t')[0]
        rxnEq = line.split('\t')[1]

        if rxnEq.count('n'):
            continue

        reaction = Reaction(rxnId)

        rxnEqRs = rxnEq.split('<=>')[0].strip()
        rxnEqPs = rxnEq.split('<=>')[1].strip()

        for comps in [rxnEqRs,rxnEqPs]:
            
            if comps.count('+'):
                comps = comps.split('+')
                for c in comps:
                    reaction = addMeToRxn(c,reaction,rxnEqPs,KeggIdModelId,idName)
            else:
                reaction = addMeToRxn(comps,reaction,rxnEqPs,KeggIdModelId,idName)
        
        # print 'Metabolites',len(reaction.Metabolites)
        reactions.append(reaction)

    universal.add_reactions(reactions)   
    print 'len(universal.reactions)',len(universal.reactions)

    return universal

def addMeToRxn(met,rxn,rxnEqPs,KeggIdModelId,idName):
    
    met = met.strip()
    num = 1
    if met.count(' '):
        num = met.split(' ')[0]
        met = met.split(' ')[1]

    if met == rxnEqPs:
        num = -num

    if met in KeggIdModelId.keys():
        ModelId = KeggIdModelId[met]
        met = Metabolite(ModelId,name = idName[met],compartment = 'c')
    else:
        met = Metabolite(met,name = idName[met],compartment = 'c')

    met_num = dict()
    met_num[met] = int(num)
    
    rxn.add_metabolites(met_num)

    return rxn

def tmp(file,coefficient = False):

    rxnIdCom = dict()
    allCom = list()

    for line in file:
        rxnId = line.split('\t')[0]
        Comps = [i.strip() for i in line.split('\t')[1].replace('<=>','+').split('+')]
        
        if coefficient :
            _Comps = copy.deepcopy(Comps)
            for i,j in enumerate(Comps):
                if j.count(' '):
                # if j.count(' ') and j.split(' ',1)[0].isalnum() :
                    j = j.split(' ')[1]
                    _Comps[i] = j
            Comps = _Comps

        rxnIdCom[rxnId] = Comps
        allCom.extend(Comps)
    allCom = list(set(allCom))

    return (rxnIdCom,allCom)

def getKeggIdModelId():

    '''In model iJO1366,metabolites's id is a series abbreviation of compoud name
    here we add the kegg reaction and metabolite needed to be identify by model shoud 
    be change the kegg id to model id '''

    f = json.load(open('./modules/iJO1366/BIGGMetInfo(1805)_203(1).json'))
    
    KeggIdModelId = dict()
    
    for i in f:
        KEGGID_list = i["database_links"]["KEGGID"] 
        modelId = i["id"]

        for _id in KEGGID_list:
            # if _id in KeggIdModelId.keys() and KeggIdModelId[_id] != modelId:
                # print _id,modelId
            KeggIdModelId[_id] = modelId

    with open('./modules/iJO1366/KeggIdModelId.json','w') as fn:
        json.dump(KeggIdModelId,fn,indent = 2)
    
    return KeggIdModelId
        

def getRxnIdEqName():

    keggRxnData = open('./keggDatebase/EnzymaticReactions.keg')
    f = open('./keggDatebase/RxnIdEqName.txt','w')

    for line in keggRxnData:
        if line.startswith('E        '):
            RxnIdEqName = line.split('E        ')[1].strip()
            RxnId = RxnIdEqName.split('  ')[0]
            EqName = RxnIdEqName.split('  ')[1]
            f.write(RxnId+'\t'+EqName+'\n')
    f.close()


def findBlock():

    '''find the block reaction by cobrapy's function-findblockreation'''

    model = cobra.io.load_json_model('./modules/iJO1366/iJO1366.json')
    model_extend = cobra.io.load_json_model('./modules/iJO1366/iJO1366_extend_rev.json')

    block = cobra.flux_analysis.find_blocked_reactions(model) #878
    block_extend = cobra.flux_analysis.find_blocked_reactions(model_extend) # 5011

    sb = [ b for b in block if b not in block_extend]
    unsb =[ b for b in block if b not in sb]
    
    print 'len(block)',len(block)
    print 'len(block_extend )',len(block_extend)

    print sb
    print 'len(sb):',sb #304
    print '-'*50

    print unsb
    print 'len(unsb)',unsb #574
    print '-'*50


def findGaps(modelFile):

    '''find gap mets by coefficient'''

    model = json.load(open(modelFile))
    # model_extend = json.load(open('./modules/iJO1366/iJO1366_extend_rev.json'))

    metsFlux = dict() #{A:{A_c:[-1,+1],A_p:[_1,]..}}
    
    for reaction in model['reactions']:
        metabolites = reaction['metabolites']
        for submet,flux in metabolites.items():

            met = submet.rsplit('_',1)[0].strip()

            if met not in metsFlux.keys():
                metsFlux[met] = dict()
                metsFlux[met][submet] = list()
                metsFlux[met][submet].append(flux)
            else:
                if submet not in metsFlux[met]:
                    metsFlux[met][submet] = list()
                metsFlux[met][submet].append(flux)


    deadEndsMets = copy.deepcopy(metsFlux)
    noProMets = dict()
    noComMets = dict()

    for met,submetFlux in metsFlux.items():
        fluxlist = reduce(lambda x,y : x + y ,submetFlux.values())

        if max(fluxlist) > 0 and min(fluxlist) < 0 :
            deadEndsMets.pop(met)
        elif max(fluxlist) < 0 :
            noProMets[met] = submetFlux
        elif min(fluxlist) > 0 :
            noComMets[met] = submetFlux


    # print 'len(metsFlux)',len(metsFlux)
    # print 'len(deadEndsMets)',len(deadEndsMets)
    # print 'len(noProMets)',len(noProMets)
    # print 'len(noComMets)',len(noComMets)
    # print noProMets

    return metsFlux

def getNosyn(modelFile):
   
    #calulate the mets's flux,get the one no flux
    model = cobra.io.load_json_model(modelFile)
    Metabolites = model.metabolites
    
    noSyn = list()
    
    f = open('./modules/iJO1366/iJO1366_extend_rev_optimizeValue.txt',"a")
    f.write('\n')
    f1 = open('./modules/iJO1366/iJO1366_extend_rev_noSyn.txt',"a")
    f1.write('\n')

    for mets in Metabolites[1116:]:
        if mets.id.startswith('add') or mets.id.startswith('No_id'):
            continue
        else:
            _model = cobra.io.load_json_model(modelFile)
            rxn = Reaction('testReaction') 
            rxn.add_metabolites({mets: -1 })
            _model.add_reaction(rxn)
            _model.change_objective('testReaction')
            val = _model.optimize().f

            f.write(mets.id+'\t'+str(val)+'\n')
            f.flush()

            if  -0.001 < val < 0.001 :
                noSyn.append(mets.id)
                f1.write(mets.id+'\t'+str(val)+'\n')
                f1.flush()


    f.close()
    f1.close()
    
    with open('./modules/iJO1366/iJO1366_extend_rev_noSyn.json','w') as fn:
        json.dump(noSyn,fn,indent = 2)
    
    print 'len(noSyn)',len(noSyn)


def index2Id(modelFile,deadEndsFile,savefile):

    '''transform the detectedDeadEnds mets's indices to ids'''
    detectDeadends = list()

    deadEnds = open(deadEndsFile)
    model = cobra.io.load_json_model(modelFile)
    Metabolites = model.metabolites

    for line in deadEnds:
        index = int(line.strip())-1
        detectDeadends.append(Metabolites[index].id)

    with open(savefile,'w') as fk:
        json.dump(detectDeadends,fk,indent = 2)
    
    return detectDeadends

def txt2json(filePath,jsonPath,ty='lis'):
    '''ty = [list,dict]'''
    f = open(filePath)
    
    if ty == 'lis':
        _list = list()
        for line in f:
            ele = line.strip().replace('\'','').strip()
            _list.append(ele)
        savefile = _list
    elif ty == 'dic':
        _dict = dict()
        seg = str(raw_input('input segmentation:'))
        for line in f:
            key = line.strip().split(seg)[0]
            val = line.strip().split(seg)[1]
            _dict[key] = val
        savefile = _dict
    else:
        print'filetype not been identified'

    with open(jsonPath,'w') as fn:
        json.dump(savefile,fn,indent = 2)
    
    return savefile

def deadEndsCla(modelFile,Deadends,comsufile,produfile):

    model = cobra.io.load_json_model(modelFile)
    metsFlux = findGaps(modelFile)

    comsu = list()
    produ = list()
    addmets = list()
    
    for submet in Deadends:

        if submet.startswith('add') or submet.startswith('No_id'):
            continue
        else:

            met = submet.rsplit('_',1)[0]
            flux = list(set(metsFlux[met][submet]))[0]
            if flux > 0 :
                produ.append(submet)
            elif flux < 0:
                comsu.append(submet)
            else:
                print submet

    print 'len(comsu)',len(comsu)
    print 'lem(produ)',len(produ)
    print 'all: ', len(comsu)+len(produ)
    
    with open(comsufile,'w') as f1:
        json.dump(comsu,f1,indent = 2)

    with open(produfile,'w') as f2:
        json.dump(produ,f2,indent = 2)

def allNoSyn(lists,filePath):

    All =  reduce(lambda x,y:x+y,lists)

    print 'before dedup:',len(All)
    print 'after dedup:',len(list(set(All)))

    dedupAll = list(set(All))
    
    with open(filePath,'w') as f:
        json.dump(dedupAll,f,indent = 2)
    
    return dedupAll

def getSmile(filePath,savepath):
    
    f = json.load(open('./modules/iJO1366/BIGGMetInfo(1805)_203(1).json'))
    ChEBIIdSmi = json.load(open('./ChEBIId_Smile.json'))
    modeId2ChEBI = dict()
    
    for met in f:
        modelId = met["id"] # is a string
        ChEBIId = met["database_links"]["CHEBI"] #is a list
        modeId2ChEBI[modelId] = ChEBIId

    metsSmi = dict()
    metsNoSmi = json.load(open(filePath))
    for mets in metsNoSmi:
        try:
            for chid in modeId2ChEBI[mets]:
                try:
                    smi = ChEBIIdSmi[chid]
                    break
                except:
                    print 'ChEBI id  %s have no smiles' % chid
            if smi:
                metsSmi[mets] = smi
            else:
                'this mets have ChEBI id but no smiles'
        except:
            print 'this mets have no ChEBI id'
    with open(savepath,'w') as fn:
        json.dump(metsSmi,fn,indent = 2)
    
    return metsSmi

def dedup(f):

    '''remove the keys has duplicated values of dict'''
    
    _f = dict()
    for i,j in f.items():
        if j not in _f.values():
            _f[i] = j
    print 'before deduplication',len(f)
    print 'After  deduplication',len(_f)
    
    return _f

def predictRxn(metsSmi):
   
    metsSmi_dedup = dedup(metsSmi)

    for mets,smi in metsSmi_dedup.items():
        if smi.count('.'): 
            continue
        BioReactor(smi,queryPsSmi=False,Ec=False,CALSIMI=True,Draw=False)

def standerdSmi(metsSmi):

    '''standerlize the smiles after NeutraliseCharges and input and output by rdkit'''
    _metsSmi = copy.deepcopy(metsSmi)

    metsList = metsSmi.keys()
    smiList = list()

    for mets,smi in metsSmi.items():
        try:
            smi = NeutraliseCharges(mtsmi(mfsmi(smi)))[0]
        except:
            smi = NeutraliseCharges(mtsmi(mfsma(smi)))[0]

        _metsSmi[mets] = smi

        if smi.count('.'):
            smi = [i.strip() for i in smi.split('.')]
            smiList += smi
        else:
           smiList.append(smi) 

    smiList = list(set(smiList))

    return _metsSmi,metsList,smiList

def validPre(metsAllFile,metsCanSynFile,metsAllNoSynFile,resultDir,Choice = 2):
    
    A = standerdSmi(dedup(json.load(open(metsAllFile))))[2]
    B = standerdSmi(dedup(json.load(open(metsCanSynFile))))[2]
    C = standerdSmi(dedup(json.load(open(metsAllNoSynFile))))[2]

    choice = {1:'atLeastOnePsInMetsAllNoSyn',
              2:'atLeastOnePsInMetsAllNoSynAndRemainedPsInMetsAll '}

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
               
                try:
                    r = NeutraliseCharges(mtsmi(mfsmi(r)))[0]
                except:
                    r = NeutraliseCharges(mtsmi(mfsma(r)))[0]

                if r not in B:
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
                        if p in C:
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
                        if p in C:
                            can = True
                        elif p in B:
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

def validPredSum(metsAllFile,metsCanSynFile,metsAllNoSynFile,resultDir,Choice = 2):

    A = standerdSmi(dedup(json.load(open(metsAllFile))))[2]
    B = standerdSmi(dedup(json.load(open(metsCanSynFile))))[2]
    C = standerdSmi(dedup(json.load(open(metsAllNoSynFile))))[2]

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
                
                if p in C:
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

    with open(os.path.join(resultDir,'preNoPro_querySmirks_%s.json' % Choice),'w') as fn:
        json.dump(preNoPro_querySmirks,fn,indent = 2)

    return preNoPro_querySmirks

def getModelSmiMets(modelMetsAll):
    
    modelMetsAll = standerdSmi(dedup(modelMetsAll))[0]
    modelSmiMets = dict()
    for mets,smi in modelMetsAll.items():
        if smi.count('.'):
            smi = [i.strip() for i in smi.split('.')]
            tmpDic = dict.fromkeys(smi,mets)
            modelSmiMets.update(tmpDic)
        else:
            modelSmiMets[smi] = mets
    return modelSmiMets
   

def metsSmi2modelId(smile,modelMetsAll):

    modelSmiMets = getModelSmiMets(modelMetsAll)

    try:
        smile = NeutraliseCharges(mtsmi(mfsmi(smile)))[0]
    except:
        smile = NeutraliseCharges(mtsmi(mfsma(smile)))[0]

    s_ModlId = modelSmiMets[smile]

    return s_ModlId

def rxnSmirks2modelId(smirks,modelMetsAll):

    rs = [i.strip() for i in smirks.split('>>')[0].split('.')]
    ps = [j.strip() for j in smirks.split('>>')[1].split('.')]

    modelSmiMets = getModelSmiMets(modelMetsAll)

    _rs = dict()
    _ps = dict()
    
    index = {'0':_rs,'1':_ps,'2':rs,'3':ps}
   
    rs_dedup = list(set(rs))
    ps_dedup = list(set(ps))

    for i,it in enumerate([rs_dedup,ps_dedup]):
        for s in it:
            s_num = index[str(i+2)].count(s)

            if i == 1:
                s_num = -s_num
            # print 'sssssss:::::',s
            try:
                s = NeutraliseCharges(mtsmi(mfsmi(s)))[0]
            except:
                s = NeutraliseCharges(mtsmi(mfsma(s)))[0]

            s_ModlId = modelSmiMets[s]
            
            index[str(i)][s_ModlId] = s_num

    _rs.update(_ps)

    return _rs

def addRxn2Model(predictRxnFile,modelFile = False):

    predictRxn = json.load(open(predictRxnFile))
    modelMetsAll = json.load(open('./modules/iJO1366/iJO1366_metsAll_smi.json'))
    
    _predictRxn = dict()

    for ps,quitem in predictRxn.items():
        
        ps_modelId = metsSmi2modelId(ps,modelMetsAll)
        
        querSmirks= dict()
        for quer,smitem in quitem.items():
            quer_modelId = metsSmi2modelId(quer,modelMetsAll)

            indexSmirks = dict()
            for i,smirks in smitem.items():
                rxn = rxnSmirks2modelId(smirks,modelMetsAll)
                indexSmirks[i] = rxn
            
            querSmirks[quer_modelId] = indexSmirks
    
        _predictRxn[ps_modelId] = querSmirks

    with open('./modules/iJO1366/iJO1366_predictRxn_modelId.json','w') as fn:
        json.dump(_predictRxn,fn,indent = 2)
    
    return _predictRxn

def gapfill(predictRxn,modelFile,saveFile):

    f = open(saveFile,'w')

    _model = cobra.io.load_json_model(modelFile)
    mets = _model.metabolites

    for ps,queryItems in predictRxn.items():
        f.write('ps'+':'+ps+'\n\n')

        for query,indexSmirks in queryItems.items():
            f.write('query'+':'+query+'\n\n')

            for i,rxndict in indexSmirks.items():
                
                model = cobra.io.load_json_model(modelFile)
                
                addrxn = dict()
                for s,val in rxndict.items():
                    addrxn[mets.get_by_id(s)] = val
                
                add_predicRxn = Reaction('add_predicRxn') 
                add_predicRxn.add_metabolites(addrxn)
                model.add_reaction(add_predicRxn)
                
                rxn = Reaction('testReaction') 
                rxn.add_metabolites({mets.get_by_id(ps): -1 })
                model.add_reaction(rxn)
                model.change_objective('testReaction')
                value = model.optimize().f

                f.write('reaction index'+str(i)+'\t'+'flux'+':'+str(value)+'\n\n')
                f.flush()
        f.write('$$$$'+'\n\n')
    f.close()


def main():
    modelDir = './modules/iJO1366/'
    modelname = 'iJO1366' 
    # modelname = 'iJO1366_extend_rev'
    modelFile = os.path.join(modelDir,'%s.json' % modelname)
    # predictRxn = json.load(open('./modules/iJO1366/iJO1366_predictRxn_modelId.json'))
    # saveFile = './modules/iJO1366/iJO1366_addpredicRxn_flux.json'
    # predictRxnFile = './preNoPro_querySmirks_2.json'

    # addRxn2Model(predictRxnFile)
    # gapfill(predictRxn,modelFile,saveFile)
    
    # f = json.load(open('./modules/iJO1366/iJO1366_predictRxn_modelId.json'))
    # _107 = f.keys()
    # _54 = list()

    # f1 = open('./modules/iJO1366/iJO1366_addpredicRxn_flux.json-bk')
    # for line in f1:
    #     if line.startswith('ps:'):
    #         c = line.split('ps:')[1].strip()
    #         _54.append(c)
    # f2 = json.load(open('./modules/iJO1366/iJO1366_noSyn+comsu+rooGaps.json'))
    
    # _107_54 = [str(i) for i in _107 if i not in _54]
   
    # _error = [j for j in _107_54 if j not in f2]
    # _right = [k for k in _107_54 if k in f2]

    # print '_107',len(_107)
    # print '_54',len(_54)
    # print '_107_54',len(_107_54)

    # print'_error',len(_error)
    # print '_right',len(_right)


    # for met,smi in C.items():
    #     if smi.count(f):
    #         print met


    # print len(allsmirk)
    # print len(list(set(allsmirk)))
    # print '-'*50
    #         allsmirk.append(smirks)
    # print len(allsmirk)


   
    # metsAllFile = os.path.join(modelDir,'%s_metsAll_smi.json' % modelname)
    # metsCanSynFile = os.path.join(modelDir,'%s_metsCanSyn_smi.json' % modelname)
    # resultDir = './preResult/'
    # validPre(metsAllFile,metsCanSynFile,metsAllNoSynFile,resultDir,Choice = 2)
    # validPredSum(metsAllFile,metsCanSynFile,metsAllNoSynFile,resultDir,Choice = 2)
    

    ChEBIIdSmi = json.load(open('./ChEBIId_Smile.json'))


    getSmile(filePath,savepath)
    print  ChEBIIdSmi[]
def script():
    pass
    '''
    
    # deadEndsFile = os.path.join(modelDir,'%s_detectDeadends.txt' % modelname)
    # savefile = os.path.join(modelDir,'%s_detectDeadends.json'% modelname)

    # comsufile = os.path.join(modelDir,'%s_comsu.json' % modelname)
    # produfile = os.path.join(modelDir,'%s_produ.json' % modelname)

    # deadEnds = index2Id(modelFile,deadEndsFile,savefile)
    # deadEndsCla(modelFile,deadEnds,comsufile,produfile)
    
    filePath = os.path.join(modelDir,'%s_noSyn+comsu+rooGaps.json' % modelname)

    model = cobra.io.load_json_model(modelFile)

    metsAll = list()
    for met in model.metabolites:

        if  met.id.startswith('add') or met.id.startswith('No_id'):
            continue
        else:
             metsAll.append(met.id)
   
    with open(os.path.join(modelDir,'%s_metsAll.json' % modelname),'w') as f:
        json.dump(metsAll,f,indent = 2)
    noSyn = json.load(open(filePath))

    metsCanSyn = [i for i in metsAll if i not in noSyn]

    with open(os.path.join(modelDir,'%s_metsCanSyn' % modelname),'w') as f1:
        json.dump(metsCanSyn,f1,indent = 2)

    print 'metsAll:',len(metsAll)

    print 'noSyn:',len(noSyn)
    print len(metsAll)-len(noSyn)
    print 'len(metsCanSyn)',len(metsCanSyn)
'''


if __name__ == '__main__':
    main()
