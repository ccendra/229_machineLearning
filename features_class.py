import functions as f
import numpy as np


class Features(object):
    """ Generator of features for one single training example. Feature is 1D vector with 136 rows. Features objects
    have the following properties:
    
    Attributes:
        name: a string representing the chemical formula of the material
        data: data frame containing all chemical information
        chemicals list: list of chemical elements
        composition: 118x1 numpy array representing chemical elemental composition
        ele: str list containing elements present in material
        stoi: float list containing stoichiometries for each element
        average atomic mass: 
    """

    def __init__(self, name, chem_list, data):
        self.name = name
        self.chem_list = chem_list
        self.data = data
        self.composition = np.zeros((len(self.chem_list), 1))
        self.weights = []
        self.ele = []
        self.stoi = []
        self.avg_mass = 0
        self.avg_column = 0
        self.avg_row = 0
        self.max_diff_an = 0
        self.avg_an = 0
        self.max_diff_ar = 0
        self.avg_ar = 0
        self.max_diff_en = 0
        self.avg_en = 0
        self.avg_s = 0
        self.avg_p = 0
        self.avg_d = 0
        self.avg_f = 0
        self.frac_s = 0
        self.frac_p = 0
        self.frac_d = 0
        self.frac_f = 0
        
    def elemental_composition(self):
        self.name = f.check_artifacts_name(self.name, self.chem_list)
        self.ele, self.stoi = f.parse_chemical_name(self.name, self.chem_list)

        for i in range(len(self.chem_list)):
            for j in range(len(self.ele)):
                if self.chem_list[i] == self.ele[j]:
                    self.composition[i, 0] = self.stoi[j]

        self.composition = self.composition/(np.sum(self.composition))
        self.weights = self.stoi/(np.sum(self.stoi))
        
        return self.composition, self.ele, self.stoi
    
    def get_rest_features(self):
        atomic_mass = self.data['atomicMass']
        atomic_mass = atomic_mass.apply(lambda x: x.replace('(', '')).apply(lambda x: x.replace(')', '')).\
            apply(lambda x: x.replace('[', '')).apply(lambda x: x.replace(']', ''))
        atomic_mass = atomic_mass.astype(float).fillna(0.0)
        
        atomic_columns = self.data['column']
        
        atomic_rows = self.data['row']
        
        atomic_number = self.data['atomicNumber']
        atomic_numbers_list = []
        
        atomic_radius = self.data['atomicRadius']
        atomic_radii_list = []
        
        atomic_electroneg = self.data['electronegativity']
        atomic_en_list = []
        
        s_electrons = self.data['s'].fillna(0.0)
        p_electrons = self.data['p'].fillna(0.0)
        d_electrons = self.data['d'].fillna(0.0)
        f_electrons = self.data['f'].fillna(0.0)
        
        for i in range(len(self.ele)):
            element = self.ele[i]
            self.avg_mass += atomic_mass.loc[element]*self.weights[i]
            self.avg_column += atomic_columns.loc[element]*self.weights[i]
            self.avg_row += atomic_rows.loc[element]*self.weights[i]
            atomic_numbers_list.append(atomic_number[element])
            self.avg_an += atomic_number.loc[element]*self.weights[i]
            atomic_radii_list.append(atomic_radius[element])
            self.avg_ar += atomic_radius.loc[element]*self.weights[i]
            atomic_en_list.append(atomic_electroneg[element])
            self.avg_en += atomic_electroneg.loc[element]*self.weights[i]
            self.avg_s += s_electrons.loc[element]*self.weights[i]
            self.avg_p += p_electrons.loc[element]*self.weights[i]
            self.avg_d += d_electrons.loc[element]*self.weights[i]
            self.avg_f += f_electrons.loc[element]*self.weights[i]
            
            den = self.avg_s + self.avg_p + self.avg_d + self.avg_f
            if den > 0:
                self.frac_s = self.avg_s/den
                self.frac_p = self.avg_p/den
                self.frac_d = self.avg_d/den
                self.frac_f = self.avg_f/den
            
        ##
        self.max_diff_an = max(atomic_numbers_list) - min(atomic_numbers_list)
        self.max_diff_ar = max(atomic_radii_list) - min(atomic_radii_list)
        self.max_diff_en = max(atomic_en_list) - min(atomic_en_list)
        
        output = np.zeros(17)
        output[0] = self.avg_mass
        output[1] = self.avg_column
        output[2] = self.avg_row
        output[3] = self.max_diff_an
        output[4] = self.avg_an
        output[5] = self.max_diff_ar
        output[6] = self.avg_ar
        output[7] = self.max_diff_en
        output[8] = self.avg_en
        output[9] = self.avg_s
        output[10] = self.avg_p
        output[11] = self.avg_d
        output[12] = self.avg_f
        output[13] = self.frac_s
        output[14] = self.frac_p
        output[15] = self.frac_d
        output[16] = self.frac_f
        
        return output