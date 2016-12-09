# !/usr/bin/python
#  _*_ coding:utf-8 _*_
# 2016/11/23
# Author:LingWu
# Email:wu_l@tib.cas.cn

import os
import copy
import json
import cobra
import rdkit
from rxnpool20161121 import *
from rdkit import Chem
from rdkit.Chem import AllChem
from cobra import Model, Reaction, Metabolite
from rdkit.Chem.AllChem import MolFromSmiles as mfsmi
from rdkit.Chem.AllChem import MolFromSmarts as mfsma
from rdkit.Chem.AllChem import MolToSmiles as mtsmi
from rdkit.Chem.AllChem import MolToSmarts as mtsma