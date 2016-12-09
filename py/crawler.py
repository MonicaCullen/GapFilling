import os
import json
from pyquery import PyQuery as pq 
from rdkit.Chem import AllChem 

def crawlerKEGGCP(keggId):

	link = 'http://www.genome.jp/dbget-bin/www_bget?-f+m+compound+%s'  % keggId
	moltxt = pq(link)

	return moltxt.text()

def iJO1366_CPSmi():

	'extract the mol of cp with keggId from kegg website ,and transfer mol to smi'
	f =json.load(open('../modules/iJO1366/BIGGMetInfo(1805)_203(1).json'))
	fk = open('../modules/iJO1366/iJO1366_CPSmiFromKEGG-1.txt','w')

	modelidSmiFromkegg = dict()
	for cp in f:
		modelid = cp["id"]
		try:
			keggId = cp["database_links"]["KEGGID"]

			for _id in keggId:
				moltxt = crawlerKEGGCP(_id)
				f1 = open('./tmp.mol' ,'w')
				f1.write('\n\n\n'+' ')
				f1.write(moltxt)
				f1.flush()
				f1.close()

				mol = AllChem.MolFromMolFile('./tmp.mol')
				smi = AllChem.MolToSmiles(mol)
				modelidSmiFromkegg[modelid] = smi
				fk.write(modelid+'\t'+smi+'\n')
				fk.flush()


		except:
			# print '{} have no keggId'.format(modelid)
			pass
	with open('../modules/iJO1366/iJO1366_CPSmiFromKEGG-1.json','w') as fn:
		json.dump(modelidSmiFromkegg,fn,indent = 2)
	fk.close()

def main():
	iJO1366_CPSmi()



if __name__ == '__main__':
	main()

