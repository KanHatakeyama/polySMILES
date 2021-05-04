"""
Class to calculate phys props by group contribution method

original code from
https://github.com/curieshicy/JRgui


"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import Counter
from operator import itemgetter
import numpy as np
import math

class JRWrapper:
    def __init__(self,temperature=298):
        self.temperature=temperature
    def compute_phys_properties(self,_SMILES):

        def return_non_duplicate_index(tuples): ##Given a list of sets return index of non_duplicate items

        ##step 1, create a new tuple, named "new_tuples"
            new_tuples = [] ##the elements are the sets
            for i in tuples:
                for j in i:
                    new_tuples.append(set(j))

        ##step 2, create a dictionary storing one to one relationship between new_tuple and old_tuple
            values = []
            for index, item in enumerate(tuples):
                if len(item) == 1:
                    values.append(index)
                else:
                    for i in [index]*len(item):
                        values.append(i)
                
            keys = [i for i in range(len(new_tuples))]
            dict_tuples = {}  ## {0:0, 1:1, 2:2, 3:3, 4:3, 5:3, 6:4, 7:4, 8:4, 9:5, 10:6, 11:7, 12:8}
            for i, j in zip(keys, values):
                dict_tuples[i] = j

        ##step 3, remove duplicates in sets terminology    
            remove_index = []
            for index_1, item in enumerate(new_tuples): ##starting from beginning
                for index_2 in range(index_1 + 1, len(new_tuples)): ##loop over the rest items
                    if len(item & new_tuples[index_2]) != 0:
                        if len(item)>len(new_tuples[index_2]):
                            remove_index.append(index_2) #indefoirx
                        elif len(item)<len(new_tuples[index_2]):
                            remove_index.append(index_1) #index
                        elif len(item)==len(new_tuples[index_2]):
                            remove_index.append(index_2) #index

            remain_sets = set(range(len(new_tuples))).difference(set(remove_index))

        ##step 4, spit out final index and length
            index_1 = [] ## [0,1,2,3,3,3,4,4]
            index_length = []

            for i in remain_sets:
                index_1.append(dict_tuples[i])
            
            ##count

            counts = Counter(index_1) ##this is a dictionary return Counter({3:3, 4:2, 0:1, 1:1, 2:1}) ##index:length

            list_counts = counts.most_common() ## convert to a list [(3,3), (4,2), (0,1), (1,1), (2,1)]

            for i in range(len(list_counts)):
                index_length.append([list_counts[i][0], list_counts[i][1]])


            index_length = sorted(index_length, key = itemgetter(0))


            return index_length

        def search_func_groups(smile): ##this is to search functional groups and print out them with numbers

            smarts = ["[$([CX2H0](=*)=*)]", "[$([CX2H1]#[!#7])]", "[$([CX2H0]#[!#7])]", "[OX2H]-[C]=O", "[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]",
          "[$([#6X3H0](=[OX1]));!$([#6X3](=[#8X1])~[#8X2]);R]=O", "[CH;D2;$(C-!@C)](=O)", "[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]", "[O;H1;$(O-!@c)]", "[#8X2H0;R;!$([#8X2H0]~[#6]=[#8])]",
          "[$([CX3H0](=[OX1]));!$([CX3](=[OX1])-[OX2]);!R]=O", "[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]", "[$([#7X3,#7X3+][!#8])](=[O])~[O-]", "[OX1H0;!$([OX1H0]~[#6X3]);!$([OX1H0]~[#7X3]~[#8])]", "[#7X2H0;R]",
          "[#7X3H1;R]", "[#7X2H1]","[#7X2H0;!R]","[#6X2]#[#7X1H0]","[NX3H2]",
          "[NX3H1;!R]", "[#7X3H0;!$([#7](~O)~O)]","[SX2H]","[#16X2H0;!R]","[#16X2H0;R]",
          "[R;CX3H1,cX3H1]", "[$([R;#6X3H0]);!$([R;#6X3H0]=[#8])]","[R;CX4H2]","[R;CX4H]","[R;CX4H0]",
          "[CX3H2]", "[!R;CX3H1;!$([CX3H1](=O))]","[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]","[CX4H3]","[!R;CX4H2]",
          "[!R;CX4H]","[!R;CX4H0]","[F]","[Cl]","[Br]",
          "[I]"]    


            tuples = []
            index_list = []
            final_index_and_length = []
            m = Chem.MolFromSmiles(str(smile))    
            for index, smart in enumerate(smarts):
                if m.HasSubstructMatch(Chem.MolFromSmarts(smart)) == True: 
                    tuples.append(m.GetSubstructMatches(Chem.MolFromSmarts(smart))) ## this is atom position
                    index_list.append(index)
           
            temp = return_non_duplicate_index(tuples) # [[0, 1], [1, 1], [3, 1], [4, 7], [5, 6], [6, 1], [7, 1], [8, 1], [9, 1]]

            for i in temp:
                final_index_and_length.append([index_list[i[0]], i[1]])

            return final_index_and_length      

        
        ##[[], [], ...[]] in total 41 nested list inside a list
        DB = [[0.0026,	0.0028,	36,	26.15,	17.78,	142.14,	136.70,	2.74E+1,	-5.57E-2,	1.01E-4,	-5.02E-8,	4.720,	2.661,	"NA",	"NA"],
                    [0.0027,	-0.0008,46,	9.20,	-11.18,	79.30,	77.71,	2.45E+1,	-2.71E-2,	1.11E-4,	-6.78E-8,	2.322,	1.155,	"NA",	"NA"],
                    [0.0020,	0.0016,	37,	27.38,	64.32,	115.51,	109.82,	7.87,	         2.01E-2,	-8.33E-6,	1.39E-9,	4.151,	3.302,	"NA",	"NA"],
                    [0.0791,	0.0077,	89,	169.09,	155.50,	-426.72,-387.87,2.41E+1,	4.27E-2,	8.04E-5,       -6.87E-8,	11.051,	19.537,	1317.23,-2.578],
                    [0.0481,	0.0005,	82,	81.10,	53.60,	-337.92,-301.95,2.45E+1,	4.02E-2,	4.02E-5,	-4.52E-8,	6.959,	9.633,	483.88,	-0.966],
                    [0.0284,	0.0028,	55,	94.97,	75.97,	-164.50,-126.27,3.04E+1,	-8.29E-2,	2.36E-4,       -1.31E-7,	"NA",	6.645,	"NA",	"NA"],
                    [0.0379,	0.0030,	82,	72.24,	36.90,	-162.03,-143.48,3.09E+1,	-3.36E-2,	1.60E-4,       -9.88E-8,	3.197,	9.093,	740.92,	-1.713],
                    [0.0741,	0.0112,	28,	92.88,	44.45,	-208.04,-189.20,2.57E+1,	-6.91E-2,	1.77E-4,	-9.88E-8,	2.406,	16.826,	2173.72,-5.057],
                    [0.0240,	0.0184,	-25,	76.34,	82.83,	-221.65,-197.37,-2.81,	         1.11E-1,	-1.16E-4,	4.94E-8,	4.490,	12.499,	3018.17,-7.314],
                    [0.0098,	0.0048,	13,	31.22,	23.05,	-138.16,-98.22,	1.22E+1,	-1.26E-2,	6.03E-5,	-3.86E-8,	5.879,	4.682,	440.24,	-0.953],
                    [0.0380,	0.0031,	62,	76.75,	61.20,	-133.22,-120.50,6.45,	         6.70E-2,	-3.57E-5,	2.86E-9,	4.189,	8.972,	340.35,	-0.350],
                    [0.0168,	0.0015,	18,	22.42,	22.23,	-132.22,-105.00,2.55E+1,	-6.32E-2,	1.11E-4,       -5.48E-8,	1.188,	2.410,	122.09,	-0.386],
                    [0.0437,	0.0064,	91,	152.54,	127.24,	-66.57,	-16.83,	2.59E+1,	-3.74E-3,	1.29E-4,	-8.88E-8,	9.679,	16.738,	"NA",	"NA"],
                    [0.0143,	0.0101,	36,	-10.50,	2.08,	-247.61,-250.83,6.82,	        1.96E-2,	1.27E-5,	-1.78E-8,	3.624,	5.909,	675.24,	-1.340],
                    [0.0085,	0.0076,	34,	57.55,	68.40,	55.52,	79.93,	8.83,	      -3.84E-3,	         4.35E-5,	-2.60E-8,	3.649,	6.528,	"NA",	"NA"],
                    [0.0130,	0.0114,	29,	52.82,	101.51,	31.65,75.61,1.18E+1,	        -2.30E-2,	1.07E-4,	-6.28E-8,	7.490,	6.930,	"NA",	"NA"],
                    ["NA",  	"NA",	"NA",	83.08,	68.91,	93.70,	119.66,	5.69,	       -4.12E-3,	1.28E-4,	-8.88E-8,	"NA",	12.169,	"NA",	"NA"],
                    [0.0255,	-0.0099,"NA",	74.60,	"NA",	23.61,	"NA",     "NA",	         "NA",	         "NA",	          "NA", 	"NA",	3.335,	"NA",	"NA"],
                    [0.0496,	-0.0101,91,	125.66,	59.89,	88.43,	89.22,	3.65E+1,	-7.33E-2,	1.84E-4,	-1.03E-7,	2.414,	12.851,	"NA",	"NA"],
                    [0.0243,	0.0109,	38,	73.23,	66.89,	-2.02,	14.07,2.69E+1,	        -4.12E-2,	1.64E-4,        -9.76E-8,	3.515,	10.788,	"NA",	"NA"],
                    [0.0295,	0.0077,	35,	50.17,	52.66,	53.47,	89.39,-1.21,	        7.62E-2,	-4.86E-5,	1.05E-8,	5.009,	6.436,	"NA",	"NA"],
                    [0.0169,	0.0074,	9,	11.74,	48.84,	123.34,	163.16,-3.11E+1,	2.27E-1,	-3.20E-4,	1.46E-7,	4.703,	1.896,	"NA",	"NA"],
                    [0.0031,	0.0084,	63,	63.56,	20.09,	-17.33,	-22.99,	3.53E+1,	-7.58E-2,	1.85E-4,	-1.03E-7,	2.360,	6.884,	"NA",	"NA"],
                    [0.0119,	0.0049,	54,	68.78,	34.40,	41.87,	33.12,	1.96E+1,	-5.61E-3,	4.02E-5,	-2.76E-8,	4.130,	6.817,	"NA",	"NA"],
                    [0.0019,	0.0051,	38,	52.10,	79.93,	39.10,	27.76,	1.67E+1,	4.81E-3,	2.77E-5,	-2.11E-8,	1.557,	5.984,	"NA",	"NA"],
                    [0.0082,	0.0011,	41,	26.73,	8.13,	2.09,	11.30,	-2.14,	        5.74E-2,	-1.64E-6,	-1.59E-8,	1.101,	2.544,	259.65,	-0.702],
                    [0.0143,	0.0008,	32,	31.01,	37.02,	46.43,	54.05,	-8.25,	        1.01E-1,	-1.42E-4,	6.78E-8,	2.394,	3.059,	-245.74,0.912],
                    [0.0100,	0.0025,	48,	27.15,	7.75,	-26.80,	-3.68,	-6.03,	         8.54E-2,	-8.00E-6,	-1.80E-8,	0.490,	2.398,	307.53,	-0.798],
                    [0.0122,	0.0004,	38,	21.78,	19.88,	8.67,	40.99,	-2.05E+1,	1.62E-1,	-1.60E-4,	6.24E-8,	3.243,	1.942,	-394.29,1.251],
                    [0.0042,	0.0061,	27,	21.32,	60.15,	79.72,	87.88,	-9.09E+1,	5.57E-1,	-9.00E-4,	4.69E-7,	-1.373,	0.644,	"NA",	"NA"],
                    [0.0113,	-0.0028,56,	18.18,	-4.32,	-9.630,	3.77,	2.36E+1,	-3.81E-2,	1.72E-4,	-1.03E-7,	-0.473,	1.724,	495.01,	-1.539],
                    [0.0129,	-0.0006,46,	24.96,	8.73,	37.97,	48.53,	-8.00,	        1.05E-1,	-9.63E-5,	3.56E-8,	2.691,	2.205,	82.28,	-0.242],
                    [0.0117,	0.0011,	38,	24.14,	11.14,	83.99,	92.36,	-2.81E+1,	2.08E-1,	-3.06E-4,	1.46E-7,	3.063,	2.138,	"NA",	"NA"],
                    [0.0141,	-0.0012,65,	23.58,	-5.10,	-76.45,	-43.96,	1.95E+1,	-8.08E-3,	1.53E-4,	-9.67E-8,	0.908,	2.373,	548.29,	-1.719],
                    [0.0189,	0.0000,	56,	22.88,	11.27,	-20.64,	8.42,	-9.09E-1,	9.50E-2,	-5.44E-5,	1.19E-8,	2.590,	2.226,	94.16,	-0.199],
                    [0.0164,	0.0020,	41,	21.74,	12.64,	29.89,	58.36,	-2.30E+1,	2.04E-1,	-2.65E-4,	1.20E-7,	0.749,	1.691,	-322.15,1.187],
                    [0.0067,	0.0043,	27,	18.25,	46.43,	82.23,	116.02,	-6.62E+1,	4.27E-1,	-6.41E-4,	3.01E-7,	-1.460,	0.636,	-573.56,2.307],
                    [0.0111,	-0.0057,27,	-0.03,	-15.78,	-251.92,-247.19,2.65E+1,	-9.13E-2,	1.91E-4,	-1.03E-7,	1.398,	-0.670,	"NA",	"NA"],
                    [0.0105,	-0.0049,58,	38.13,	13.55,	-71.55,-64.31,	3.33E+1,	-9.63E-2,	1.87E-4,	-9.96E-8,	2.515,	4.532,	625.45,	-1.814],
                    [0.0133,	0.0057,	71,	66.86,	43.43,	-29.48,	-38.06,	2.86E+1,	-6.49E-2,	1.36E-4,	-7.45E-8,	3.603,	6.582,	738.91,	-2.038],
                    [0.0068,	-0.0034,97,	93.84,	41.69,	21.06,	5.74,	3.21E+1,	-6.41E-2,	1.26E-4,	-6.87E-8,	2.724,	9.520,	809.55,	-2.224]] 

        self.NoA = Chem.AddHs(Chem.MolFromSmiles(str(_SMILES))).GetNumAtoms()
        self.MW = Descriptors.MolWt(Chem.AddHs(Chem.MolFromSmiles(str(_SMILES))))
        self.LogP = Descriptors.MolLogP(Chem.AddHs(Chem.MolFromSmiles(str(_SMILES))))
        self.MR = Descriptors.MolMR(Chem.AddHs(Chem.MolFromSmiles(str(_SMILES))))

        double_lists = search_func_groups(_SMILES)


        entry_index_by_users = []
        entry_data_by_users = []

        for item in double_lists:
            entry_index_by_users.append(item[0])
            entry_data_by_users.append(item[1])
            
        fiveteen_columns = [] ##length  = 15*len(entry_index_by_users)
        for index, data in zip(entry_index_by_users, entry_data_by_users):
            for i in range(15):
                if DB[index][i] == np.nan:
                    temp = np.nan
                else:
                    temp = data*DB[index][i]
                fiveteen_columns.append(temp)

        Tc = []
        Pc = []
        Vc = []
        Tb = []
        Tm = []
        Hfor = []
        Gf = []
        Cpa = []
        Cpb = []
        Cpc = []
        Cpd = []
        Hfus = []
        Hvap = []
        Ya = []
        Yb =[]        
        fc = fiveteen_columns ## short hand
        for i in range(len(entry_index_by_users)):
            Tc.append(fc[i*15])
            Pc.append(fc[i*15 + 1])
            Vc.append(fc[i*15 + 2])
            Tb.append(fc[i*15 + 3])
            Tm.append(fc[i*15 + 4])
            Hfor.append(fc[i*15 + 5])
            Gf.append(fc[i*15 + 6])
            Cpa.append(fc[i*15 + 7])
            Cpb.append(fc[i*15 + 8])
            Cpc.append(fc[i*15 + 9])
            Cpd.append(fc[i*15 + 10])
            Hfus.append(fc[i*15 + 11])
            Hvap.append(fc[i*15 + 12])
            Ya.append(fc[i*15 + 13])
            Yb.append(fc[i*15 + 14])
        try:
            self._BoilingPoint = 198.2 + sum(Tb)
        except:
            self._BoilingPoint = np.nan        
        try:
            self._MeltingPoint = 122.5 + sum(Tm)
        except:
            self._MeltingPoint = np.nan
        try:
            self._CriticalTemp =  (sum(Tb) + 198.2)/(0.584 + 0.965*sum(Tc) - sum(Tc)**2)
        except:
            self._CriticalTemp = np.nan
        
        try:
            self._CriticalPress = 1./(0.113 + 0.0032*float(self.NoA) - sum(Pc))**2
        except:
            self._CriticalPress = np.nan
        try:
            self._CriticalVolume = 17.5 + sum(Vc)
        except:
            self._CriticalVolume = np.nan
        try:
            self._EnthalpyForm = 68.29 + sum(Hfor)
        except:
            self._EnthalpyForm = np.nan
        try:
            self._GibbsEnergy = 53.88 + sum(Gf)
        except:
            self._GibbsEnergy = np.nan
        try:
            self._HeatCapacity = (sum(Cpa) - 37.93) + (sum(Cpb) + 0.210)*float(self.temperature) + (sum(Cpc) - 3.91*10**(-4))*float(self.temperature)**2 + (sum(Cpd) + 2.06*10**(-7))*float(self.temperature)**3 
        except:
            self._HeatCapacity = np.nan
        try:
            self._EnthalpyVap = 15.30 + sum(Hvap)
        except:
            self._EnthalpyVap = np.nan
        try:
            self._EnthalpyFus = -0.88 + sum(Hfus)
        except:
            self._EnthalpyFus = np.nan
        try:
            self._LiquidVisco = float(self.MW)*math.exp((sum(Ya) - 597.82)/float(self.temperature) + sum(Yb) - 11.202)
        except:
            self._LiquidVisco = np.nan

        try:
            self._CrystalSolub_1 = 10**(0.8 - float(self.LogP) - 0.01*(sum(Tm)+122.5 - 273.15 - 25.))*1000.*float(self.MW)
        except:
            self._CrystalSolub_1 = np.nan
        try:
            self._CrystalSolub_2 = 10**(0.5 - float(self.LogP) - 0.01*(sum(Tm)+122.5 - 273.15 - 25.))*1000.*float(self.MW)
        except:
            self._CrystalSolub_2 = np.nan
        try:
            self._AmorphSolub_1 = 10**(0.8 - float(self.LogP) - 0.01*(sum(Tm)+122.5 - 273.15 - 25.)) *1000.*float(self.MW)*math.exp((sum(Hfus)-0.88)*(sum(Tm) + 122.5 - float(self.temperature))*float(self.temperature)/(sum(Tm) + 122.5)**2/(2.479*float(self.temperature)/298.))
        except:
            self._AmorphSolub_1 = np.nan
        try:
            self._AmorphSolub_2 = 10**(0.5 - float(self.LogP) - 0.01*(sum(Tm)+122.5 - 273.15 - 25.)) *1000.*float(self.MW)*math.exp((sum(Hfus)-0.88)*(sum(Tm) + 122.5 - float(self.temperature))*float(self.temperature)/(sum(Tm) + 122.5)**2/(2.479*float(self.temperature)/298.))
        except:
            self._AmorphSolub_2 = np.nan
            

        labels=["BoilingPoint", "MeltingPoint", "CriticalTemp", "CriticalPress", "CriticalVolume",
                "EnthalpyForm", "GibbsEnergy", "HeatCapacity", "EnthalpyVap", "EnthalpyFus",
                "LiquidVisco", "CrystalSolub_1", "CrystalSolub_2", "AmorphSolub_1", "AmorphSolub_2"]
        labels=["JR_"+i for i in labels]
        vals=[self._BoilingPoint, self._MeltingPoint, self._CriticalTemp, self._CriticalPress, self._CriticalVolume,
                self._EnthalpyForm, self._GibbsEnergy, self._HeatCapacity, self._EnthalpyVap, self._EnthalpyFus,
                self._LiquidVisco, self._CrystalSolub_1, self._CrystalSolub_2, self._AmorphSolub_1, self._AmorphSolub_2] 
        return {k:v for k,v in zip(labels,vals)}