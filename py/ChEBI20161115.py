import os
import json

comps = open('../chebiData/ChEBI_complete.sdf').read().split('$$$$')

def getIdSmile():

    '''get the compound smiles ,generate a dict chebi Id as the key and smiles as the value'''

    ChEBIIdSmile = dict()

    for comp in comps:
        
        _idlist = list()
        compAll = comp.split('\n\n')
        smile = False
        
        for item in compAll:
            if item.count('> <ChEBI ID>'):
                _id = item.split('> <ChEBI ID>')[1].strip().split('CHEBI:')[1]
                _idlist.append(_id)
            elif item.count('> <Secondary ChEBI ID>'):
                subid = item.split('> <Secondary ChEBI ID>')[1].strip()
                subid = [i.strip() for i in subid.split('CHEBI:') if i]
                _idlist += subid
            elif item.count('> <SMILES>'):
                smile = item.split('> <SMILES>')[1].strip()
            else:
                pass
            if smile and _idlist:
                tmpdic = dict.fromkeys(_idlist,smile)
                ChEBIIdSmile.update(tmpdic)

    print len(ChEBIIdSmile)
    with open('../chebiData/ChEBIIdSmile.json','w') as fn:

        json.dump(ChEBIIdSmile,fn,indent = 2)

    return ChEBIIdSmile

def getSmiName():
    
    '''get the compound name from CHEBI SDF file,generate a dict -smile as the key and name as the value'''

    ChEBISmiName = dict()
    
    for comp in comps:

        smiName = dict()

        compAll = comp.split('\n\n')

        for item in compAll:
            if item.count('> <SMILES>'):
                smile = item.split('> <SMILES>')[1].strip()
            elif item.count('> <ChEBI Name>'):
                name = item.split('> <ChEBI Name>')[1].strip()
                smiName['ChEBI Name'] = name
            elif item.count('> <IUPAC Names>'):
                IUPAC_Names = item.split('> <IUPAC Names>')[1].strip()
                smiName['IUPAC_Names'] = IUPAC_Names
            elif item.count('> <Synonyms>'):
                Synonyms = item.split('> <Synonyms>')[1].strip().split('\n')
                smiName['Synonyms'] = Synonyms
            else:
                pass
        if smile and smiName:
            ChEBISmiName[smile] = smiName
    
    with open('../chebiData/ChEBISmiName.json','w') as fn:
        json.dump(ChEBISmiName,fn,indent = 2)

    return ChEBISmiName

def chebismi_NameId():
    
    ChEBISmiNameId = dict()

    for comp in comps:

        smiNameId = dict()
        _idlist = list()
        smile = False

        compAll = comp.split('\n\n')

        for item in compAll:

            if item.count('> <SMILES>'):
                smile = item.split('> <SMILES>')[1].strip()
            elif item.count('> <ChEBI ID>'):
                _id = item.split('> <ChEBI ID>')[1].strip().split('CHEBI:')[1]
                _idlist.append(_id)
            elif item.count('> <Secondary ChEBI ID>'):
                subid = item.split('> <Secondary ChEBI ID>')[1].strip()
                subid = [i.strip() for i in subid.split('CHEBI:') if i]
                _idlist += subid
            elif item.count('> <ChEBI Name>'):
                name = item.split('> <ChEBI Name>')[1].strip()
                smiNameId['ChEBI Name'] = name
            elif item.count('> <IUPAC Names>'):
                IUPAC_Names = item.split('> <IUPAC Names>')[1].strip()
                smiNameId['IUPAC_Names'] = IUPAC_Names
            elif item.count('> <Synonyms>'):
                Synonyms = item.split('> <Synonyms>')[1].strip().split('\n')
                smiNameId['Synonyms'] = Synonyms
            else:
                pass
            if _idlist:
                smiNameId['ChEBI ID'] = _idlist
            if smile :
                ChEBISmiNameId[smile] = smiNameId
    with open('../chebiData/ChEBISmi_NameId.json','w') as fn:
        json.dump(ChEBISmiNameId,fn,indent=2)

def main():
    # getIdSmile()
    # chebismi_NameId()
    f = json.load(open('../chebiData/ChEBISmi_NameId.json'))
    print len(f.keys())

    



if __name__ == '__main__':
    main()