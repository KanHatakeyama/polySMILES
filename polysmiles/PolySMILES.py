from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
from rdkit import Chem
import numpy as np

from read_smiles import read_smiles
from polymer_graph_helper import rename_polymer_nodes,update_QQ_nodes,update_Q_nodes,fragmentate_units,draw_chem_graph
from polymer_graph_helper import graph_to_dict,process_molecular_weight


descriptor_names = [descriptor_name[0] for descriptor_name in Descriptors._descList]
descriptor_calculation = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
default_calculator=descriptor_calculation.CalcDescriptors


class PolySMILES:

    def __init__(self,
                calculator=default_calculator,
                descriptor_names=descriptor_names,
                cap_atom="H",
                fill_nan=0
                ):
        self.calculator=calculator
        self.descriptor_names=descriptor_names
        self.cap_atom=cap_atom
        self.fill_nan=fill_nan

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
        mol = Chem.MolFromSmiles(smiles)
        descriptor=np.array(self.calculator(mol))

        if self.fill_nan is not None:
            descriptor=np.nan_to_num(descriptor,nan=self.fill_nan)

        return descriptor

    
    def smiles_to_weighted_descriptors(self,
                                       smiles,
                                       default_n=50
                                      ):

        pol_dict=self.smiles_to_dict(smiles,calculate_descriptor=True)        

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
        desc_dict={"total MW":total_mw}
        temp_dict={k:v for k,v in zip(self.descriptor_names,average_descriptor)}
        desc_dict.update(temp_dict)

        return desc_dict