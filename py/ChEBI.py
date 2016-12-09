import os
import json

def getSmile():

    comps = open('./ChEBI_complete.sdf').read().split('$$$$')

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
            else:
                pass

            if smile and _idlist:
                tmpdic = dict.fromkeys(_idlist,smile)
                ChEBIIdSmile.update(tmpdic)

    print len(ChEBIIdSmile)
    with open('./ChEBIId_Smile.json','w') as fn:

        json.dump(ChEBIIdSmile,fn,indent = 2)

    return ChEBIIdSmile

    
def main():
    getSmiles()

if __name__ == '__main__':
    main()