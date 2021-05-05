from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import networkx as nx
import re
import copy

from .write_smiles import write_smiles


#***** graph processing funcs ***********
SHARE_PARAMS=["mw","mn","n","pdi","d"]
DUMMY_ATOM="Y"
DUMMY_ATOM_WEIGHT=rdMolDescriptors._CalcMolWt(Chem.MolFromSmiles("["+DUMMY_ATOM+"]"))
MAX_ATOMS=10000  


def draw_chem_graph(g):
    """
    draw chemicals from networkX object
    
    Parameters
    ----------
    g : networkX object

    """
    pos = nx.spring_layout(g)  
    node_labels = nx.get_node_attributes(g,'polymer')
    nx.draw_networkx_labels(g, pos, labels = node_labels,font_size=10, font_color="b")
    node_labels = nx.get_node_attributes(g,'element')
    nx.draw_networkx_labels(g, pos, labels = node_labels,font_size=16, font_color="r")
    nx.draw_networkx_labels(g, pos, font_size=9, font_color="r")
    nx.draw_networkx_nodes(g, pos, node_size=100, node_color="w")
    nx.draw_networkx_edges(g, pos, width=1)
    
    
def get_element(g,i):
    return g.nodes[i]["element"]

def parse_polymer_label(txt):
    txt=txt.replace("Pol","").replace("Q","").replace(" ","").replace(":","").replace("pdi","d")

    prop_dict={}
    for i in txt.split(","):
        if i !="":
            prop,value=i.split("=")
            prop_dict[prop]=value

    return prop_dict


def rename_polymer_nodes(g):
    for i in g.nodes:
        label= get_element(g,i)
        # find polymer-related nodes
        if re.findall('Q|Pol', label):
            # set parsed parameters to "polymer" attribute
            g.nodes[i]["polymer"]=parse_polymer_label(label)
            # reset node name
            g.nodes[i]["element"]="Q"

            

# split Q-Q nodes
def update_QQ_nodes(g):

    cut_edge_list=[]

    #explore "Q" nodes
    for i in nx.get_node_attributes(g,'polymer').keys():
        #find neighboring "Q" nodes
        for j in list(g.neighbors(i)):
            label=get_element(g,j)
            if label!="Q":
                continue        

            #at this step, i and j nodes are both "Q" nodes
            #write connection information
            for m,n in ((i,j),(j,i)):
                if "connect" in g.nodes[m]["polymer"].keys():
                     g.nodes[n]["polymer"]["connect"]=g.nodes[m]["polymer"]["connect"]


            #delete connection
            g.remove_edge(i, j)
            cut_edge_list.append((i,j))
            
    return cut_edge_list            


#share molecular weight info in the graphs


def update_Q_nodes(g):
    #explore "Q" nodes
    for i in nx.get_node_attributes(g,'polymer').keys():
        #search "Q" nodes connected in the same graph
        for j in range(i):
            if i==j or get_element(g,j)!="Q":
                continue

            #check connection
            if nx.has_path(g, source=i, target=j):

                #add nw and mn info
                for param in SHARE_PARAMS: 
                    for m,n in ((i,j),(j,i)):
                        if param in g.nodes[m]["polymer"].keys():
                             g.nodes[n]["polymer"][param]=g.nodes[m]["polymer"][param]


        #add implicit connection information (if no info, this indicates "block" connection)
        if "connect" not in g.nodes[i]["polymer"].keys():
            g.nodes[i]["polymer"]["connect"]="block"
            
            
def fragmentate_units(g):
    #search Q nodes
    for i in nx.get_node_attributes(g,'polymer').keys():

        neighbor_node_ids=list(g.neighbors(i))
        if len(neighbor_node_ids)==1:
            continue
        elif len(neighbor_node_ids)>2:
            raise ValueError("error! Q nodes must be connected to <= 2 nodes.")

        #duplicate a Q node
        q_node1=i
        q_node2=(q_node1)+MAX_ATOMS
        
        g.add_node(q_node2)

        #set attributes
        for k,v in g.nodes[q_node1].items():
            g.nodes[q_node2][k]=v

        #reconnect nodes
        g.remove_edge(q_node1, neighbor_node_ids[0])
        g.add_edge(q_node2,neighbor_node_ids[0])
        #g.edges[q_node2,neighbor_node_ids[0]]["order"]=1
        
        

def graph_to_dict(g):

    #extract info from subgraphs
    total_dict={}
    for graph_num,sub_g in enumerate([g.subgraph(c) for c in nx.connected_components(g)]):
        #sub_g=copy.deepcopy(sub_g)
        #extract polymer info

        #explore Q nodes
        pol_dict={}
        q_node_list=nx.get_node_attributes(sub_g,'polymer').keys()
        for num,i in enumerate(q_node_list):

            #if the unit is polymeric
            if len(q_node_list)>1:
                if num==0:
                    # copy polymer info (just use the 1st Q node because info except for "connect" is the same in a graph)
                    pol_dict=copy.copy(sub_g.nodes[i]["polymer"])
                else:
                    #just update connect info.
                    pol_dict["connect"]=pol_dict["connect"]+"-"+sub_g.nodes[i]["polymer"]["connect"]

            #extact smiles info

            #sub_g.nodes[i]["element"]="Y"
            #sub_g.nodes[i].pop("polymer")


        #try:
        pol_dict["SMILES"]=write_smiles(sub_g)      
        #except:
        #    pass

        total_dict[graph_num]=pol_dict

    return total_dict



#*****funcs to process dict data************

def calc_molecular_weight(sm):
    sm=sm.replace("Q",DUMMY_ATOM)
    mol = Chem.MolFromSmiles(sm)
    mw=rdMolDescriptors._CalcMolWt(mol)
    mw=mw-DUMMY_ATOM_WEIGHT*sm.count(DUMMY_ATOM)
    return mw

def process_molecular_weight(pol_dict):
    #calculate molecular weight for each unit
    for unit in pol_dict.keys():
        sm=pol_dict[unit]["SMILES"]
        mw=calc_molecular_weight(sm)
        pol_dict[unit]["mw_unit"]=mw

        if sm.count("Q")<2:
            pol_dict[unit]["type"]="monomeric"
        else:
            pol_dict[unit]["type"]="polymeric"

        #set molecular weights for monomeric units    
        if pol_dict[unit]["type"]=="monomeric":
            pol_dict[unit]["mn"]=mw
            pol_dict[unit]["mw"]=mw     


        #calculate Mw, Mn ,etc for polymeric units
        if pol_dict[unit]["type"]=="polymeric":

            #calculate Mn from n
            if "n" in pol_dict[unit].keys():
                try:
                    n=float(pol_dict[unit]["n"])
                except:
                    raise ValueError("error! n is not number")

                pol_dict[unit]["mn"]=float(pol_dict[unit]["mw_unit"])*n

            #calculate Mn from Mw and d
            if ("mw" in pol_dict[unit].keys()) and ("d" in pol_dict[unit].keys()):
                try:
                    mw=float(pol_dict[unit]["mw"])
                    d=float(pol_dict[unit]["d"])
                except:
                    print("error! mw or d are not number")      

                pol_dict[unit]["mn"]=mw/d

            #calculate Mw from Mn and d
            if ("mn" in pol_dict[unit].keys()) and ("d" in pol_dict[unit].keys()):
                try:
                    mn=float(pol_dict[unit]["mn"])
                    d=float(pol_dict[unit]["d"])
                except:
                    print("error! mn or d are not number")      

                pol_dict[unit]["mw"]=mn*d

            #calculate d from Mw and Mn
            if ("mn" in pol_dict[unit].keys()) and ("mw" in pol_dict[unit].keys()):
                try:
                    mn=float(pol_dict[unit]["mn"])
                    mw=float(pol_dict[unit]["mw"])
                except:
                    print("error! mn or mw are not number")      

                pol_dict[unit]["d"]=mw/mn

            #calculate n from Mn
            if "mn" in pol_dict[unit].keys():
                try:
                    mn=float(pol_dict[unit]["mn"])
                except:
                    raise ValueError("error! mn is not number")

                pol_dict[unit]["n"]=mn/float(pol_dict[unit]["mw_unit"])

