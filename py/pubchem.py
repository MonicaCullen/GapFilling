import os
import json

def pubSmi_NameId():

	pubchemFile = open('../pubchemData/Compound_000175001_000200000.sdf').read()
	comps = pubchemFile.split('$$$$')
	
	pubSmiNameId = dict()

	for comp in comps:
		
		comAll = comp.split('\n\n')
		NameId = dict()
		CAN_SMILES = False
		
		for item in comAll:
			if item.count('> <PUBCHEM_COMPOUND_CID>'):
				comId = item.split('> <PUBCHEM_COMPOUND_CID>')[1].strip()
				NameId['pubChemId'] = comId
			elif item.count('> <PUBCHEM_IUPAC_NAME>'):
				IUPAC_NAME = item.split('> <PUBCHEM_IUPAC_NAME>')[1].strip()
				NameId['IUPAC_NAME'] = IUPAC_NAME
			elif item.count('> <PUBCHEM_OPENEYE_CAN_SMILES>'):
				CAN_SMILES = item.split('> <PUBCHEM_OPENEYE_CAN_SMILES>')[1].strip()
			else:
				pass
			if CAN_SMILES:
				pubSmiNameId[CAN_SMILES] = NameId
	print  len(pubSmiNameId.keys())
	with open('../pubchemData/pubSmi_NameId.json','w') as fn:
		json.dump(pubSmiNameId,fn,indent = 2)
	
	return pubSmiNameId


def main():
	pubSmi_NameId()

if __name__ == '__main__':
	main()