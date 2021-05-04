from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
from rdkit import Chem
import numpy as np
import pandas as pd


from read_smiles import read_smiles
from polymer_graph_helper import rename_polymer_nodes,update_QQ_nodes,update_Q_nodes,fragmentate_units,draw_chem_graph
from polymer_graph_helper import graph_to_dict,process_molecular_weight

from AutoDescriptor import AutoDescriptor



class PolySMILES:

    def __init__(self,
                calculator=AutoDescriptor(),
                cap_atom="H",
                 
                ):
        self.calculator=calculator
        self.cap_atom=cap_atom
        self.dict_mode=True
        
    def auto(self,smiles_list):
        if type(smiles_list) is type(""):
            smiles_list=[smiles_list]
        res_dict={}
        for i,smiles in enumerate(smiles_list):
            try:
            #if True:
                res_dict[i]=self.smiles_to_weighted_descriptors(smiles)
            except:
                print(i, " error!",  smiles)
        return pd.DataFrame.from_dict(res_dict).T


    def smiles_to_dict(self,
                       smiles,
                       graph_mode=False,
                       calculate_descriptor=False
                      ):
        g=read_smiles(smiles)

        #process polymer info described in Q nodes
        rename_polymer_nodes(g)

        #cut polymeric units connected with Q-Q nodes
        cut_edge_list=update_QQ_nodes(g)

        if graph_mode:
            return g,cut_edge_list
        
        #share info with Q nodes connected in the same polymer units
        update_Q_nodes(g)

        #fragmentate units
        fragmentate_units(g)

        #to dict
        pol_dict=graph_to_dict(g)
        
        #calcualte molecular weight
        process_molecular_weight(pol_dict)
        
        self.g=g
        self.pol_dict=pol_dict
        
        #calculate descriptor
        if calculate_descriptor:
            for unit in pol_dict.keys():
                sm=pol_dict[unit]["SMILES"]
                pol_dict[unit]["descriptor"]=self._calc_descriptor(sm)
                
        return pol_dict
    
    
    def _calc_descriptor(self,smiles):
        smiles=smiles.replace("Q",self.cap_atom)
        res_df=self.calculator(smiles)
        res_dict={k:v for k,v in zip(res_df.columns,res_df.values[0])}
        res_dict.pop("SMILES")
        
        if self.dict_mode:
            return res_dict
        else:
            return list(res_dict.values())
    
    def smiles_to_weighted_descriptors(self,
                                       smiles,
                                       default_n=50
                                      ):

        temp_mode=self.dict_mode
        self.dict_mode=False
        pol_dict=self.smiles_to_dict(smiles,calculate_descriptor=True)        
        self.dict_mode=temp_mode
        
        mw_array=[]
        #calculate molecular weight of each unit
        for unit in pol_dict.keys():

            # 1) use actual Mw
            if "mw" in pol_dict[unit].keys():
                mw_array.append(pol_dict[unit]["mw"])
            # 2) or use actual Mn
            elif "mn" in pol_dict[unit].keys():
                mw_array.append(pol_dict[unit]["mn"])
            # 3) is there is no info, just estimate Mw assuming n=50
            else:
                mw_array.append(pol_dict[unit]["mw_unit"]*default_n)

        #calculate weight ratio
        mw_array=np.array(mw_array)
        total_mw=np.sum(mw_array)
        mw_ratio=mw_array/total_mw

        #calculate descriptors
        desc_array=np.array([pol_dict[unit]["descriptor"] for unit in pol_dict.keys()])
        average_descriptor=np.dot(desc_array.T,mw_ratio)

        #return as dict
        desc_dict={"SMILES":smiles,"total MW":total_mw}
        temp_dict={k:v for k,v in zip(self.calculator.descriptor_names,average_descriptor)}
        desc_dict.update(temp_dict)

        return desc_dict