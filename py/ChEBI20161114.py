import os
import json

comps = open('./chebiData/ChEBI_complete.sdf').read().split('$$$$')

def getIdSmile():

    ChEBIIdSmile = dict()

    for comp in comps:
        
        _idlist = list()
        compAll = comp.split('\n\n')
        
        (_id,subid,smile) = (False,False,False)
        for item in compAll:
            if item.count('> <ChEBI ID>'):
                _id = item.split('> <ChEBI ID>')[1].strip().split('CHEBI:')[1]
                _idlist.append(_id)
            elif item.count('> <Secondary ChEBI ID>'):
                subid = item.split('> <Secondary ChEBI ID>')[1].strip()
                subid = [i.strip() for i in subid.split('CHEBI:') if i]
                _idlist += subid
            elif item.count('> <SMILES>'):
                smile = True
                smile = item.split('> <SMILES>')[1].strip()
            elif item.count('> <ChEBI Name>'):
                name = item.split('> <ChEBI Name>')[1].strip()
            # elif item.count('> <Definition>'):
            #     Definition = item.split('> <Definition>')[1].strip()
            # elif item.count('> <InChIKey>'):
            #     InChIKey = item.split('> <InChIKey>')[1].strip()
            # elif item.count('> <InChI>'):
            #     InChI = item.split('> <InChI>')[1].strip()
            # elif item.count('> <Formulae>'):
            #     Formulae = item.split('> <Formulae>')[1].strip()
            # elif item.count('> <Charge>'):
            #     Charge = item.split('> <Charge>')[1].strip()
            # elif item.count('> <Mass>'):
            #     Mass = item.split('> <Mass>')[1].strip()
            # elif item.count('> <Monoisotopic Mass>'):
            #     Monoisotopic_Mass = item.split('> <Monoisotopic Mass>')[1].strip()
            # elif item.count('> <IUPAC Names>'):
            #     IUPAC_Names = item.split('> <IUPAC Names>')[1].strip()
            elif item.count('> <Synonyms>'):
                Synonyms = item.split('> <Synonyms>')[1].strip().split('\n')
            else:
                pass

            if smile and _idlist:
                tmpdic = dict.fromkeys(_idlist,smile)
                ChEBIIdSmile.update(tmpdic)

    print len(ChEBIIdSmile)
    with open('./ChEBIId_Smile.json','w') as fn:

        json.dump(ChEBIIdSmile,fn,indent = 2)

    return ChEBIIdSmile

def getSmiName():
    

    ChEBISmiName = dict()
    
    for comp in comps:

        smiName = dict()

        compAll = comp.split('\n\n')

        (_id,subid,smile) = (False,False,False)

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
        with open('./chebiData/ChEBISmiName.json','w') as fn:
            json.dump(ChEBISmiName,fn,indent = 2)

def main():
    # getSmiles()
    getSmiName()

if __name__ == '__main__':
    main()