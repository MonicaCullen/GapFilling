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
from rxnpool20161206 import *
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
    
    f = json.load(open('../modules/iJO1366/BIGGMetInfo(1805)_203(1).json'))
    ChEBIIdSmi = json.load(open('../chebiData/ChEBIIdSmile.json'))
    keggMetsSmi = json.load(open('../modules/iJO1366/iJO1366_CPSmiFromKEGG.json'))
    modeId2ChEBI = dict()
    
    for met in f:
        modelId = met["id"] # is a string
        ChEBIId = met["database_links"]["CHEBI"] #is a list
        modeId2ChEBI[modelId] = ChEBIId

    metsSmi = dict()
    metsNoSmi = json.load(open(filePath))
    kegg = 0
    for mets in metsNoSmi:

        try:
            chebiHasSmi = False

            for chid in modeId2ChEBI[mets]:

                smi = False

                try:
                    smi = ChEBIIdSmi[chid]
                    chebiHasSmi = True
                except:
                    # print 'ChEBI id  %s have no smiles' % chid
                    pass
            
                if smi:
                    metsSmi[mets] = smi
                else:
                    # 'this mets have ChEBI id  %s but no smiles' % chid
                    pass

            if not chebiHasSmi:
                try:
                    smi = keggMetsSmi[mets]
                    metsSmi[mets] = smi
                    kegg += 1
                    print kegg,'%s :' % mets,smi
                except:
                    # print 'KEGG also have no smiles of this mets'
                    continue
        except:
            # print 'this mets have no ChEBI id'
            pass
    with open(savepath,'w') as fn:
        json.dump(metsSmi,fn,indent = 2)
    
    return metsSmi

def missingExchangeRxn(metsCanSynFile,metsNonSynFile):
    '''find the missing exchange reactions between metsCanSyn and metsNosyn'''
    metsCanSynDic = json.load(open(metsCanSynFile))
    metsNonSynDic = json.load(open(metsNonSynFile))

    metsCanSyn = [i.rsplit('_',1)[0].strip() for i in metsCanSynDic.keys() ]
    metsNonSyn = [i.rsplit('_',1)[0].strip() for i in metsNonSynDic.keys() ]

    fakeNonSyn = [i for i in metsNonSyn if i in metsCanSyn]
    
    canSyn_mets_metSub = dict()

    for metsub in metsCanSynDic.keys():
        mets = metsub.rsplit('_',1)[0]
        if mets not in canSyn_mets_metSub.keys():
            canSyn_mets_metSub[mets] = list()
        canSyn_mets_metSub[mets].append(metsub)
        canSyn_mets_metSub[mets] = list(set(canSyn_mets_metSub[mets]))
    
    # nonSyn_mets_metSub = dict()
    fakeNonSynAndCanSyn = dict()
    for metsub in metsNonSynDic.keys():
        mets = metsub.rsplit('_',1)[0]
       
        if mets in fakeNonSyn:
            metPre = canSyn_mets_metSub[mets]

            if mets not in fakeNonSynAndCanSyn.keys():
                fakeNonSynAndCanSyn[mets] = dict()
            
            fakeNonSynAndCanSyn[mets][metsub] = metPre

    return fakeNonSynAndCanSyn
    
def addExchangeRxn(modelFile,missingExRxn):

    # f =open('./addExchangeRxn_fluxChange_onebyone.txt','w')
    f =open('../modules/iJO1366/original_model/preResult/calFakeNoSynCpFlux.txt','w')

    fakeNonSynAndCanSyn = missingExRxn

    fakeNoSynCpAddRxns = dict()
    for mets,metsub in fakeNonSynAndCanSyn.items():

        fakeNoSynCpAddRxns[mets] = dict()
        
        model = cobra.io.load_json_model(modelFile)
        
        noSyn = metsub.keys()
        
        if len(noSyn) == 1:

            m = len(model.reactions)
           
            sub = noSyn[0]
            
            fakeNoSynCpAddRxns[mets][sub] = dict()

            '''sub's _p existed and have flux,identify there be or not be a rxn to transport sub out'''
            
            (mets_e,addRxnOut,model) = addMetsRxnOut(sub,model)

            model.change_objective(addRxnOut) #改变目标函数为代谢物的分泌反应
            
            beforeFlux = model.optimize().f #填gap之前的流量值
            fakeNoSynCpAddRxns[mets][sub]['ori_flux'] = beforeFlux

            # 添加_p-->_e 的交换反应
            mets_p = model.metabolites.get_by_id(sub.replace('_e','_p'))

            addRxnName = 'add_exchangeRxn_%s_p-e' % mets
            fakeNoSynCpAddRxns[mets][sub][addRxnName] = {'%s_p' % sub:-1,'%s__e' % sub:1}

            exchangeRxn = cobra.Reaction(addRxnName)
            
            exchangeRxn.add_metabolites({mets_p:-1,mets_e:1})
            
            model.add_reaction(exchangeRxn)

            model.change_objective(addRxnOut)
            
            afterFlux = model.optimize().f
            fakeNoSynCpAddRxns[mets][sub]['new_flux'] = afterFlux
            
            n = len(model.reactions)

            # f.write('ModelRxns:  '+ str(m) + '\n')
            # f.write('ModelRxns:  ' + str(n) + '\n')
            # f.write('add   rxn:  '+ str(n-m) + '\n')
            # f.write('noSyn  CP:  ' + sub + '\n')
            # f.write('beforFlux:  ' + str(beforeFlux) + '\n') 
            # f.write('afterFlux:  '+ str(afterFlux) + '\n') 
            # f.write('-'*50 + '\n')
            # f.flush()
            f.write('*'*70 + '\n')
            f.write(mets + ' '*(32-len(mets)) + 'ori_flux_e :' +str(beforeFlux) + '\n'*2)
            f.write('-'*70 +'\n')
            f.write(addRxnName +' '*(32-len(addRxnName))+ 'new_flux_e :'+ str(afterFlux) + '\n'*2)

            # print fakeNoSynCpAddRxns
            # print '-'*50

        else:

            m = len(model.reactions)

            for sub in noSyn:

                fakeNoSynCpAddRxns[mets][sub] = dict()
                
                if sub.endswith('_e'):

                    (mets_e,addRxnOut_e,model) = addMetsRxnOut(sub,model)

                    mets_c = model.metabolites.get_by_id(sub.replace('_e','_c')) #修改@1215（把_c 写成了 _p）

                    model.change_objective(addRxnOut_e) #改变目标函数为代谢物的分泌反应

                    beforeFlux_e = model.optimize().f #填gap之前的流量值
                    fakeNoSynCpAddRxns[mets][sub]['ori_flux'] = beforeFlux_e


                elif sub.endswith('_p'):
                    
                    (mets_p,addRxnOut_p,model) = addMetsRxnOut(sub,model)
                    
                    model.change_objective(addRxnOut_p) #改变目标函数为代谢物的分泌反应

                    beforeFlux_p = model.optimize().f #填gap之前的流量值
                    fakeNoSynCpAddRxns[mets][sub]['ori_flux'] = beforeFlux_p


                else:
                    print sub
            addRxnName_p_e = 'add_exchangeRxn_%s_p-e' % mets
            addRxnName_c_p = 'add_exchangeRxn_%s_c-p' % mets

            fakeNoSynCpAddRxns[mets][mets + '_p'][addRxnName_c_p] = {'%s_c' % mets:-1,'%s_p' % mets:1}
            fakeNoSynCpAddRxns[mets][mets + '_e'][addRxnName_c_p] = {'%s_c' % mets:-1,'%s_p' % mets:1}
            fakeNoSynCpAddRxns[mets][mets + '_e'][addRxnName_p_e] = {'%s_p' % mets:-1,'%s_e' % mets:1}

            exchangeRxn_p_e = cobra.Reaction('add_exchangeRxn_%s_p-e' % mets)
            exchangeRxn_c_p = cobra.Reaction('add_exchangeRxn_%s_c-p' % mets)

            exchangeRxn_p_e.add_metabolites({mets_p:-1,mets_e:1})
            exchangeRxn_c_p.add_metabolites({mets_c:-1,mets_p:1})

            model.add_reactions([exchangeRxn_p_e,exchangeRxn_c_p])

            model.change_objective(addRxnOut_p)
            
            afterFlux_p = model.optimize().f
            
            model.change_objective(addRxnOut_e)

            afterFlux_e = model.optimize().f

            fakeNoSynCpAddRxns[mets][mets + '_p']['new_flux'] = afterFlux_p
            fakeNoSynCpAddRxns[mets][mets + '_e']['new_flux'] = afterFlux_e


            n = len(model.reactions)

            # f.write('ModelRxns  :  ' + str(m) +'\n')
            # f.write('ModelRxns  :  ' + str(n) + '\n')
            # f.write('add     rxn:  ' + str(n-m)+ '\n')
            # f.write('noSyn    CP:  ' + str(noSyn) + '\n')
            # f.write('beforFlux_e:  ' + str(beforeFlux_e) + '\n')
            # f.write('afterFlux_e:  ' + str(afterFlux_e) + '\n')
            # f.write('beforFlux_p:  '+ str(beforeFlux_p)+ '\n')
            # f.write('afterFlux_p:  '+ str(afterFlux_p) + '\n')
            # f.write( '-'*50 + '\n')
            # f.flush()
            f.write('*'*70 + '\n')
            f.write(str(mets) + ' '*(32-len(str(mets))) + 'ori_flux_p :' +str(beforeFlux_p) + '\n'*2)
            f.write(str(mets) + ' '*(32-len(str(mets))) + 'ori_flux_e :' +str(beforeFlux_e) + '\n'*2)

            f.write('-'*70 +'\n')
            f.write(addRxnName_c_p +' '*(32-len(addRxnName_c_p))+ 'new_flux_p :'+ str(afterFlux_p) + '\n'*2)
            f.write(addRxnName_p_e +' '*(32-len(addRxnName_p_e))+ 'new_flux_e :'+ str(afterFlux_e) + '\n'*2)
            
            f.flush()

            # print fakeNoSynCpAddRxns

            # print '-'*50

    f.close
            
    with open('../modules/iJO1366/original_model/preResult/iJO1366_fakeNoSynCpAddRxns.json','w') as fn:
        json.dump(fakeNoSynCpAddRxns,fn,indent = 2)

def addMetsRxnOut(sub,model):
    
    subRxnout = '%s --> ' % sub #创建代谢物的分泌反应

    transportOut = False #默认模型中没有代谢物的分泌反应

    mets = model.metabolites.get_by_id(sub) #获取代谢物在模型中的对象

    rxns = list(mets.reactions) #获取代谢物参与的反应

    for rxn in rxns: #循环反应，查找是否已经含有分泌反应

        if rxn.reaction == subRxnout: #如果含有分泌反应，则赋值，标记找到分泌反应

            addRxnOut = rxn

            transportOut = True

    if not transportOut: #如果没有找到分泌反应，则添加分泌反应

        addRxnOut = cobra.Reaction('add_outRxn_%s' % sub)

        addRxnOut.add_metabolites({mets:-1})

        model.add_reaction(addRxnOut)

    return (mets,addRxnOut,model)

def preCanSynAndNoSyn():

    CanSyn = json.load(open('../modules/iJO1366/iJO1366_metsCanSyn_smi(kegg5).json'))
    NonSyn = json.load(open('../modules/iJO1366/iJO1366_noSyn+comsu+rooGaps_smi(kegg13).json'))
    fakeNonSyn = json.load(open('../modules/iJO1366/iJO1366_fakeNonSynAndItsCanSyn.json'))
    
    _NonSyn = copy.deepcopy(NonSyn)
    
    print 'len(CanSyn)',len(CanSyn)
    print 'len(NonSyn)',len(NonSyn)

    for metsub,smi in NonSyn.items():
        mets = metsub.rsplit('_',1)[0]
        if mets in fakeNonSyn.keys():
            _NonSyn.pop(metsub)
            CanSyn[metsub] = smi
    
    print 'len(CanSyn)',len(CanSyn)
    print 'len(_NonSyn)',len(_NonSyn)
    CanSyn_Smi_Names = dict()
    NonSyn_Smi_Names = dict()

    for mets,smi in CanSyn.items():
        if smi not in CanSyn_Smi_Names.keys():
            CanSyn_Smi_Names[smi] = list()
        CanSyn_Smi_Names[smi].append(mets)

    for mets,smi in _NonSyn.items():
        if smi not in NonSyn_Smi_Names.keys():
            NonSyn_Smi_Names[smi] = list()
        NonSyn_Smi_Names[smi].append(mets)
   
    print 'len(CanSyn_Smi_Names)',len(CanSyn_Smi_Names)
    print 'len(NonSyn_Smi_Names)',len(NonSyn_Smi_Names)

    with open('../modules/iJO1366/iJO1366_metsCanSyn+fakeNoSyn.json','w') as fn:
        json.dump(CanSyn,fn,indent = 2)
    with open('../modules/iJO1366/iJO1366_metsNonSyn-fakeNoSyn.json','w') as fn:
        json.dump(_NonSyn,fn,indent = 2)
    
    with open('../modules/iJO1366/iJO1366_metsCanSyn+fakeNoSyn_Smi-Names.json','w') as fn:
        json.dump(CanSyn_Smi_Names,fn,indent = 2)
    with open('../modules/iJO1366/iJO1366_metsNonSyn-fakeNoSyn_smi-Names.json','w') as fn:
        json.dump(NonSyn_Smi_Names,fn,indent = 2)

def CanSyn2NoSynPredict():

    CanSyn_Smi_Names = json.load(open('../modules/iJO1366/iJO1366_metsCanSyn+fakeNoSyn_Smi-Names.json'))

    for smi in CanSyn_Smi_Names.keys():
        if smi.count('.'): 
            print 'smi contain "."',smi
            continue
        BioReactor(smi,queryPsSmi=False,Ec=False,Draw=False)

def dedup(f):

    '''remove the keys has duplicated values of dict'''
    
    _f = dict()
    for i,j in f.items():
        if j not in _f.values():
            _f[i] = j
    print 'before deduplication',len(f)
    print 'After  deduplication',len(_f)
    
    return _f

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

def validPre(metsCanSyn,metsNonSyn,resultDir,Choice = 2):

 # # '''因为有部分化合物是以盐的形式存在，所以Smiles里存在'.'的情况,需要将其分开处理再去匹配，否则一些反应物
    # # 或者生成物因为smiles不能完全匹配的问题，被认为是不存在于CanSyn或者NonSyn'''
   
    # metsCanSynFile = json.load(open('../modules/iJO1366/original_model/iJO1366_metsCanSyn+fakeNoSyn_Smi-Names.json'))
    # metsNonSynFile = json.load(open('../modules/iJO1366/original_model/iJO1366_metsNonSyn-fakeNoSyn_Smi-Names.json'))
    
    # metsCanSyn = [ NeutraliseCharges(i)[0] for i in metsCanSynFile.keys()]

    # _metsCanSyn = copy.deepcopy(metsCanSyn)
    # for met in metsCanSyn:
    #     if met.count('.'):
    #         mets = met.split('.')
    #         _metsCanSyn = _metsCanSyn +mets
    
    # metsCanSyn = list(set(metsCanSyn))
    # _metsCanSyn = list(set(_metsCanSyn))

    # print 'metsCanSyn',len(metsCanSyn)
    # print '_metsCanSyn',len(_metsCanSyn)

    # metsNonSyn = [ NeutraliseCharges(i)[0] for i in metsNonSynFile.keys()]

    # _metsNonSyn = copy.deepcopy(metsNonSyn)
    # for met in metsNonSyn:
    #     if met.count('.'):
    #         mets = met.split('.')
    #         _metsNonSyn = _metsNonSyn +mets

    
    # metsNonSyn =  list(set(metsNonSyn))
    # _metsNonSyn =  list(set(_metsNonSyn))
    # _metsNonSyn.remove('NCCO')
    # _metsNonSyn.remove('CCN')

    # print 'metsNonSyn',len(metsNonSyn)
    # print '_metsNonSyn 1 ',len(_metsNonSyn)
    # #因为标准化之后。一些化合物在去掉构型之后Smi是相同的，因此会同时存在于CanSyn和NonSyn,会影响下面的Smi匹配
    # for i in _metsCanSyn:
    #     if i in _metsNonSyn:
    #         _metsNonSyn.remove(i)

    # print '_metsNonSyn 2 ',len(_metsNonSyn)

    # resultDir = '../modules/iJO1366/original_model/preResult/'
    # validPre(_metsCanSyn,_metsNonSyn,resultDir,Choice = 2)
    # preNoPro_querySmirks = validPredSum(_metsCanSyn,_metsNonSyn,resultDir,Choice = 2)
    # print len(preNoPro_querySmirks)

    
    B = metsCanSyn ; C = metsNonSyn

    choice = {1:'atLeastOnePsInMetsAllNoSyn',
              2:'atLeastOnePsInMetsAllNoSynAndRemainedPsInMetsAll '}

    f = open(os.path.join(resultDir,'%s-20161215.txt' % choice[Choice]),'w')

    for filename in os.listdir(resultDir):

        if filename.startswith('2016'):

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

                    r_stand = NeutraliseCharges(r)

                    if not r_stand:
                        r_stand = r
                    else:
                        r_stand = r_stand[0]

                    if r_stand not in B:

                        can = False
                        _resultDict.pop(str(i))
                        break
                    else:
                        pass

                if can:
                    if Choice == 1:
                        can = False
                        for p in ps:

                            p_stand = NeutraliseCharges(p)

                            if not p_stand:
                                p_stand = p
                            else:
                                p_stand = p_stand[0]
                            
                            if p_stand in C:
                                can = True
                            else:
                                pass
                        if not can:
                            _resultDict.pop(str(i))

                    elif Choice == 2:
                        remained = True
                        for p in ps:

                            p_stand = NeutraliseCharges(p)

                            if not p_stand:
                                p_stand = p
                            else:
                                p_stand = p_stand[0]
                            
                            n = 0
                            if p_stand in C:
                                n += 1
                            elif p_stand in B: #modified@20161213 p p_stand
                                pass
                            else:
                                remained = False
                        
                        if n != 1 or (not remained):
                            _resultDict.pop(str(i))
                        print 

            output = filename+'\t'+'before:'+str(len(resultDict.keys()))+'\t'+'after:'+str(len(_resultDict))+'\n'
            f.write(output)
            f.flush()

            with open(os.path.join(resultDir,filename,'result_%s-20161215.json' % Choice),'w') as fn:
                json.dump(_resultDict,fn,indent = 2)

    f.close()

def validPredSum(metsCanSyn,metsNonSyn,resultDir,Choice = 2):

    B = metsCanSyn ; C = metsNonSyn

    preNoPro_querySmirks = dict()

    for filename in os.listdir(resultDir):

        if filename.startswith('2016'):

            for name in os.listdir(os.path.join(resultDir+filename)):

                path = resultDir+filename+'/'+name
                
                if name == 'result.txt':
                    result_txt = open(path).read()
                    query = mtsmi(mfsmi(result_txt.split('\n',1)[0].split(':',1)[1].strip()))
                
                if name == 'result_%s-20161215.json' % Choice: 
                    resultDict = json.load(open(path))
        
            for i,item in resultDict.items():
                
                smirks = item["Smikrs"]
                rs = smirks.split(">>")[0].split(".")
                ps = smirks.split(">>")[1].split(".")

                for p in ps :

                    p_stand = NeutraliseCharges(p)

                    if not p_stand:
                        p_stand = p
                    else:
                        p_stand = p_stand[0]

                    if p_stand in C:
                        if p_stand not in preNoPro_querySmirks.keys():
                            preNoPro_querySmirks[p_stand] = dict()
                            preNoPro_querySmirks[p_stand][query] = dict()
                            preNoPro_querySmirks[p_stand][query][i] = smirks
                    
                        else:
                            if query in preNoPro_querySmirks[p_stand].keys():
                                preNoPro_querySmirks[p_stand][query][i] = smirks
                            else:
                                preNoPro_querySmirks[p_stand][query] = dict()
                                preNoPro_querySmirks[p_stand][query][i] = smirks
                    else:
                        pass

    with open(os.path.join(resultDir,'preNoPro_querySmirks_%s-20161215.json' % Choice),'w') as fn:
        json.dump(preNoPro_querySmirks,fn,indent = 2)

    return preNoPro_querySmirks

def preRxn2ModelRxn(Draw = False):
    preNoPro_querySmirks = json.load(open('../modules/iJO1366/original_model/preResult/preNoPro_querySmirks_2-20161215.json'))
    metsCanSyn = json.load(open('../modules/iJO1366/original_model/iJO1366_metsCanSyn+fakeNoSyn_standSmi-Names.json'))
    metsNonSyn = json.load(open('../modules/iJO1366/original_model/iJO1366_metsNonSyn-fakeNoSyn_standSmi-Names.json'))
    metsNonSyn_extend = json.load(open('../modules/iJO1366/extend_model/iJO1366_extend_rev_noSyn.json'))
    
    dotCanSynSmi = [i for i in metsCanSyn.keys() if i.count('.')]
    dotNonSynSmi = [i for i in metsNonSyn.keys() if i.count('.')]
    dotSmi = dotCanSynSmi + dotNonSynSmi
        
    dotSmi_modelId = dict()
    for Smi in dotSmi:
        try:
            Smi_modelId = metsCanSyn[Smi]
        except:
            Smi_modelId = metsNonSyn[Smi]

        Smi_list = Smi.split('.')
        for s in Smi_list:
            if s in metsCanSyn.keys():
                Smi_modelId = metsCanSyn[s]
            dotSmi_modelId[s] = Smi_modelId[0].rsplit('_',1)[0]
    # n = 0
    ps_query_index_rxnDic = dict()

    for ps,psitem in preNoPro_querySmirks.items():
                
        ps = NeutraliseCharges(ps)[0]
        ps_modelId = metsNonSyn[ps][0].rsplit('_',1)[0]

        if Draw:
            ImageDir = '../modules/iJO1366/original_model/preResult/preNoPro_Image/'
            if not os.path.exists(ImageDir):
                os.mkdir(ImageDir)

            f = open(os.path.join(ImageDir,'./tmp.smi'),'w')
            f.write(ps)
            f.flush()
            f.close()
            rxnSmiFile = os.path.join(ImageDir,'./tmp.smi')
            rxnImagefile = os.path.join(ImageDir,'%s.png' % (ps_modelId))
            DrawImage(rxnSmiFile,rxnImagefile,800,400)

        query_index_rxnDic = dict()

        for query,rxnitems in psitem.items():
            query = NeutraliseCharges(query)[0]
            query_modelId = metsCanSyn[query][0].rsplit('_',1)[0]

            index_rxnDic = dict()

            for index,smirks in rxnitems.items():
                
                rxn = dict() #一个反应创建一个字典，value为系数值，消耗为负，生成为正

                rxn_rs = smirks.split('>>')[0].split('.')
                rxn_ps = smirks.split('>>')[1].split('.')

                for r in rxn_rs:
                    r = NeutraliseCharges(r)[0]

                    if r in dotSmi_modelId.keys():
                        r_modelId = dotSmi_modelId[r]
                    else:
                        r_modelId = metsCanSyn[r][0].rsplit('_',1)[0]

                    if r_modelId not in rxn.keys():
                        rxn[r_modelId] = -1
                    else:
                        rxn[r_modelId] -= 1

                for p in rxn_ps:
                    p = NeutraliseCharges(p)[0]

                    if p in dotSmi_modelId.keys():
                        p_modelId = dotSmi_modelId[p]
                    else:
                        try :
                            p_modelId = metsCanSyn[p][0].rsplit('_',1)[0]
                        except:
                            p_modelId = metsNonSyn[p][0].rsplit('_',1)[0]

                    if p_modelId not in rxn.keys():
                        rxn[p_modelId] = 1
                    else:
                        rxn[p_modelId] += 1

                if rxn:
                    index_rxnDic[index] = rxn
                
                if Draw:
                    f = open(os.path.join(ImageDir,'./tmp.smi'),'w')
                    f.write(smirks)
                    f.flush()
                    f.close()
                    rxnSmiFile = os.path.join(ImageDir,'./tmp.smi')
                    rxnImagefile = os.path.join(ImageDir,'%s_%s_%s_rxn.png' % (ps_modelId,query_modelId,index))
                    DrawImage(rxnSmiFile,rxnImagefile,800,400)

            query_index_rxnDic[query_modelId] = index_rxnDic

        ps_query_index_rxnDic[ps_modelId] = query_index_rxnDic
    with open('../modules/iJO1366/original_model/preResult/psModelId_queryModelId_index_rxnDic-20161215.json','w') as fn:
        json.dump(ps_query_index_rxnDic,fn,indent = 2)
    print 'len(ps_query_index_rxnDic.keys())',len(ps_query_index_rxnDic.keys())

def modelRxn2txt():
   
    modelDir = '../modules/iJO1366/original_model/'
    
    filename = 'preResult/psModelId_queryModelId_index_rxnDic-20161215.json'
    file = os.path.join(modelDir,filename)
    psModelId_queryModelId_index_rxnDic = json.load(open(file))
    
    filename = 'iJO1366_metsCanSyn+fakeNoSyn.json'
    file = os.path.join(modelDir,filename)
    metsCanSyn = json.load(open(file))
    canSyn = metsCanSyn.keys()

    filename = 'iJO1366_metsNonSyn-fakeNoSyn.json'
    file = os.path.join(modelDir,filename)
    metsNonSyn = json.load(open(file))
    nonSyn = metsNonSyn.keys()

    f = open(os.path.join(modelDir,'preResult/psModelId_queryModelId_index_rxnDic-20161215.txt'),'w')
    f.write('*'*50+'\n')
    f.write('\t'*5)
    f.write('c'+'\t'+'p'+'\t'+'e'+'\n')

    for ps,qsitem in psModelId_queryModelId_index_rxnDic.items():
        f.write(ps +'(p)' +' '*(13-len(ps))+'N'+'\t')

        ps_c = ps + '_c'
        ps_p = ps + '_p'
        ps_e = ps + '_e'

        psset = [ps_c,ps_p,ps_e]
        for m,p in enumerate(psset):
            if p in nonSyn:
                f.write('1' + '\t')
            else:
                f.write('0' + '\t')

        f.write('\n'*2) 
        for qs,irxn in qsitem.items():
            f.write('-'*50+'\n'+qs +'(q)' +' '*(13-len(qs))+'Y'+'\t'+'\n'*2)

            for i,rxn in irxn.items():
                f.write('reaction '+str(i) + '\n')
                for cp,coef in rxn.items():
                    if cp == ps:
                        continue
                    else:
                        if coef < 0:
                            f.write(cp +'(r)' +' '*(13-len(cp))+'Y'+'\t')
                        else:
                            f.write(cp +'(p)' +' '*(13-len(cp))+'Y'+'\t')

                        cp_c = cp+'_c'
                        cp_p = cp+'_p'
                        cp_e = cp+'_e'

                        cpset = (cp_c,cp_p,cp_e)
                        for c in cpset:
                            if c in canSyn:
                                f.write('1' + '\t')
                            else:
                                f.write('0' + '\t')
                        f.write('\n')
                f.write('\n') 
        f.write('end'+'\n'+'*'*50+'\n'*2) 
    f.close()



def addpreModelRxn():

    modelDir = '../modules/iJO1366/original_model/'
    
    filename = 'preResult/psModelId_queryModelId_index_rxnDic-20161215.json'
    file = os.path.join(modelDir,filename)
    psModelId_queryModelId_index_rxnDic = json.load(open(file))
    
    filename = 'iJO1366_metsCanSyn+fakeNoSyn.json'
    file = os.path.join(modelDir,filename)
    metsCanSyn = json.load(open(file))
    canSyn = metsCanSyn.keys()

    filename = 'iJO1366_metsNonSyn-fakeNoSyn.json'
    file = os.path.join(modelDir,filename)
    metsNonSyn = json.load(open(file))
    nonSyn = metsNonSyn.keys()

    synCpAddRxns = dict()
    for ps,qsitem in psModelId_queryModelId_index_rxnDic.items():
        
        ps_c = ps + '_c'
       
        addRxns = dict()
        

        for qs,indexRxn in qsitem.items():

            for i,rxn in indexRxn.items():

                addRxnName = 'add_%s_%s_%s' %(ps,qs,i)

                addRxns[addRxnName] = dict()
                addRxns[addRxnName]['addmets'] = list()
                
                if ps_c not in nonSyn:
                    addRxns[addRxnName]['addmets'].append(ps_c)


                mainRxn = dict()
                coRxns = list()
                
                for cp,coef in rxn.items():
                    if cp == ps:
                        mainRxn[ps_c] = coef
                    else:
                        cp_c = cp +'_c'
                        
                        if cp_c not in canSyn:

                            addRxns[addRxnName]['addmets'].append(cp_c)
                            
                            cp_p = cp + '_p'

                            if cp_p not in canSyn:

                                addRxns[addRxnName]['addmets'].append(cp_p)

                                coRxn_e2p = dict()
                                rxnName = 'addCoRxn_%s2%s' % (cp_e,cp_p)
                                rxn = {cp_p:-1,cp_c:1}
                                coRxn_e2p[rxnName] = rxn
                                coRxns.append(coRxn_e2p)
                            
                            coRxn_p2c = dict()
                            rxnName = 'addCoRxn_%s2%s' % (cp_p,cp_c)
                            rxn = {cp_p:-1,cp_c:1}
                            coRxns.append(coRxn_p2c)
                        
                        else:
                            mainRxn[cp_c] = coef

                addRxns[addRxnName]['mainRxn'] = mainRxn
                addRxns[addRxnName]['coRxns'] = coRxns 
        
        synCpAddRxns[ps_c]=addRxns
    with open('../modules/iJO1366/original_model/iJO1366_synCpAddRxns.json','w') as fn:
        json.dump(synCpAddRxns,fn,indent = 2)
    print len(synCpAddRxns)

def calNonSynCpFlux():

    synCpAddRxns = json.load(open('../modules/iJO1366/original_model/preResult/iJO1366_synCpAddRxns.json'))
    model = cobra.io.load_json_model('../modules/iJO1366/original_model/ijO1366.json')
    
    f = open('../modules/iJO1366/original_model/preResult/calNonSynCpFlux.txt-bk','w')
   
    calNonSynCpFlux = dict()
   
    for cp,rxnitem in synCpAddRxns.items():
        f.write('*'*50+'\n')
        n =1
        
        for addrxn,item in rxnitem.items():
            _model = copy.deepcopy(model)

            #先判断添加该预测反应 是否需要添加化合物，如果需要先把化合物添加齐
            if  item["addmets"]:
                for met in item["addmets"]:
                    try:
                        _met = _model.metabolites.get_by_id(met.replace('_c','_p'))
                    except:
                        _met = _model.metabolites.get_by_id(met.replace('_c','_e'))
                    meta = cobra.Metabolite(met,formula = _met.formula,name = _met.name,charge = _met.charge,compartment = 'c')
                    _model.add_metabolites(meta)

            #添加完化合物(有些_c化合物并不存在于原模型中)，先计算添加反应前目标化合物的流量
            (cp_mets,addRxnOut,_model) = addMetsRxnOut(cp,_model)
            if n>0:
                _model.change_objective(addRxnOut)
                original_flux = _model.optimize().f
                f.write(cp + ' '*(24-len(cp)) +'ori_flux :' + str(original_flux) + '\n')
                f.write('-'*50+'\n')
                n -= 1

                calNonSynCpFlux[cp] = dict()
                calNonSynCpFlux[cp]['ori_flux'] = original_flux

            
            #添加新反应
            reaction = cobra.Reaction(addrxn)
            for met,coef in item['mainRxn'].items():
                reaction.add_metabolites({_model.metabolites.get_by_id(met):coef})

            _model.add_reaction(reaction)
            _model.change_objective(addRxnOut)
            new_flux = _model.optimize().f
            f.write(addrxn + ' '*(24-len(addrxn)) +'new_flux :' + str(new_flux) + '\n'*2)
            f.flush()
            calNonSynCpFlux[cp][addrxn] = new_flux
    f.close()
    with open('../modules/iJO1366/original_model/preResult/calNonSynCpFlux.json','w') as fn:
        json.dump(calNonSynCpFlux,fn,indent = 2)

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
    try:
        s_ModlId = modelSmiMets[smile]
        return s_ModlId
    except:
        return None

    rs = [i.strip() for i in smirks.split('>>')[0].split('.')]
    ps = [j.strip() for j in smirks.split('>>')[1].split('.')]

    modelSmiMets = getModelSmiMets(modelMetsAll)

    _rs = dict()
    _ps = dict()
    
    index = {'0':_rs,'1':_ps,'2':rs,'3':ps}
   
    rs_dedup = list(set(rs))
    ps_dedup = list(set(ps))

    for i,it in enumerate([rs_dedup,ps_dedup]):
        have = True
        for s in it:
            s_num = index[str(i+2)].count(s)

            if i == 1:
                s_num = -s_num
            # print 'sssssss:::::',s
            try:
                s = NeutraliseCharges(mtsmi(mfsmi(s)))[0]
            except:
                s = NeutraliseCharges(mtsmi(mfsma(s)))[0]

            try:
                s_ModlId = modelSmiMets[s]
                index[str(i)][s_ModlId] = s_num
            except:
                have = False
    if have:
        _rs.update(_ps)
        return _rs
    else:
        return None

def addRxn2Model(predictRxnFile,modelFile = False):

    predictRxn = json.load(open(predictRxnFile))
    modelMetsAll = json.load(open('./modules/iJO1366/iJO1366_metsAll_smi-test.json'))
    
    _predictRxn = dict()

    for ps,quitem in predictRxn.items():
        
        ps_modelId = metsSmi2modelId(ps,modelMetsAll)

        if ps_modelId:
        
            querSmirks= dict()
            for quer,smitem in quitem.items():
                
                quer_modelId = metsSmi2modelId(quer,modelMetsAll)
                
                if quer_modelId:

                    indexSmirks = dict()
                    
                    for i,smirks in smitem.items():
                       
                        rxn = rxnSmirks2modelId(smirks,modelMetsAll)

                        if rxn:
                            indexSmirks[i] = rxn
                    
                    querSmirks[quer_modelId] = indexSmirks
    
        _predictRxn[ps_modelId] = querSmirks

    with open('./modules/iJO1366/iJO1366_predictRxn_modelId-test2.json','w') as fn:
        json.dump(_predictRxn,fn,indent = 2)
    
    return _predictRxn

def gapfill(predictRxn,modelFile,saveFile):

    f = open(saveFile,'w')

    _model = cobra.io.load_json_model(modelFile)
    mets = _model.metabolites

    for ps,queryItems in predictRxn.items():
        f.write('ps'+':'+ps+'\n')

        for query,indexSmirks in queryItems.items():
            f.write('query'+':'+query+'\n')

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

                f.write('reaction index'+str(i)+'\t'+'flux'+':'+str(value)+'\n')
                f.flush()
        f.write('$$$$'+'\n\n')
    f.close()

def fluxCompare():
    filePath = open('../modules/iJO1366/iJO1366_optimizeValue.txt')
    optFlux = dict()
    for line in filePath:
        ps = line.split('\t')[0]
        flux = line.split('\t')[1].strip()
        optFlux[ps] = flux

    f = open('../modules/iJO1366/iJO1366_addpredicRxn_flux.json-test').read()
    psflux = f.split('$$$$')
    psfluxDic = dict()
    for ps in psflux:
        fluxlist = list()
        linelist = ps.split('\n')
        for line in linelist:
            if line.count('ps'):
                print ps
                print '-'*50
                p = line.split('ps:')[1].strip()
                
            elif line.count('flux'):
                flux = line.split('flux:')[1].strip()
                print flux
                fluxlist.append(flux)
        psfluxDic[p] = fluxlist

    can = list()
    for ps,fluxlist in psfluxDic.items():
        if max(fluxlist) >0.01 and max(fluxlist) > optFlux[ps]:
            can.append(ps)
    print len(can)
def writeNoSynFile():

    file = json.load(open('../modules/iJO1366/original_model/preResult/iJO1366_synCpAddRxns.json'))
    file1 = json.load(open('../modules/iJO1366/original_model/preResult/calNonSynCpFlux.json'))
    
    f = open('../modules/iJO1366/original_model/preResult/iJO1366_synCpAddRxns.txt','w')
    for cp ,addRxn in file.items():
        f.write('*'*110+'\n')
        f.write('metabolite :' + ' '*(32-len('metabolite :'))+cp +' '*(50-len(cp))+'ori_flux:'+str(file1[cp]["ori_flux"])+'\n')
        f.write('-'*110 + '\n')

        for rxnName,it in addRxn.items():
            if it["addmets"]:
                pass
                # f.write('add mets :' +' '*(32-len('add mets :'))+ it["addmets"][0] + '\n')

            rs = list()
            ps = list()

            for m,coef in it['mainRxn'].items():
                if coef < 0 :
                    if abs(coef) != 1: 
                        rs.append(str(abs(coef)) + ' ' + m)
                    else:
                        rs.append(m)
                else:
                    if abs(coef) != 1: 
                        ps.append(str(abs(coef)) + ' ' + m)
                    else:
                        ps.append(m)

            rxn = ' + '.join(rs) + ' --> '+' + '.join(ps)
            f.write(rxnName + ' '*(32-len(rxnName)) +rxn + ' '*(50-len(rxn))+'new_flux:'+str(file1[cp][rxnName])+'\n'*2)

def writeFakeNoSynFile():

    file = json.load(open('../modules/iJO1366/original_model/preResult/iJO1366_fakeNoSynCpAddRxns.json'))
    
    f = open('../modules/iJO1366/original_model/preResult/iJO1366_fakeNoSynCpAddRxns.txt','w')
    
    for met,subit in file.items():

        for cp,cpit in subit.items():

            f.write('*'*100+'\n')
            f.write('metabolite :' + ' '*(32-len('metabolite :'))+cp +' '*(35-len(cp))+'ori_flux:'+str(cpit["ori_flux"])+'\n')
            f.write('-'*100 + '\n')

            for rxnName,it in cpit.items():
                if rxnName.startswith('add'):
                    for c,coef in it.items():
                        if coef < 0:
                            rs = c
                        if coef > 0 :
                            ps = c
                    rxn = rs + ' --> '+ ps
                    f.write(rxnName + ' '*(32-len(rxnName)) +rxn + ' '*(35-len(rxn))+'new_flux:'+str(cpit['new_flux'])+'\n'*2)


    
def main():
    # addpreModelRxn()
    # calNonSynCpFlux() 
    # preRxn2ModelRxn(Draw = True)
    writeFakeNoSynFile()

if __name__ == '__main__':
    # main()
    # modelFile = '../modules/iJO1366/original_model/iJO1366.json'
    # missingExRxn = json.load(open('../modules/iJO1366/original_model/iJO1366_fakeNonSynAndItsCanSyn.json'))
    # addExchangeRxn(modelFile,missingExRxn)
    resultDir = ('../modules/iJO1366/original_model/preResult/')

    for filename in os.listdir(resultDir):

        if filename.startswith('2016'):

            for name in os.listdir(os.path.join(resultDir+filename)):

                path = resultDir+filename+'/'+name
                
                if name == 'result.txt':
                    result_txt = open(path).read()
                    query = mtsmi(mfsmi(result_txt.split('\n',1)[0].split(':',1)[1].strip()))
                    query_stand = NeutraliseCharges(query)[0]
                    if query_stand == 'OCC1OC(O)C(O)C(O)C1O':
                        print filename
                        break





'''
    metsNonSyn_extend = json.load(open('../modules/iJO1366/extend_model/iJO1366_extend_rev_noSyn.json'))
    metsNoSyn = json.load(open('../modules/iJO1366/original_model/iJO1366_metsNonSyn-fakeNoSyn.json'))
    print len(metsNoSyn.keys())
    print len(metsNonSyn_extend)

    diff = [i for i in metsNoSyn.keys() if  i in metsNonSyn_extend]
    de_diff = set(list(diff))

    diff_ = [i.rsplit('_',1)[0] for i in diff]
    de_diff_ = list(set(diff_))
    # print diff 
    print len(diff)
    print len(de_diff)

    print len(diff_)
    print len(de_diff_)

'''



