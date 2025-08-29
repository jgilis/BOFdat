"""
Protein
=======

This module generates BOFsc for the 20 amino acids contained in proteins.

"""
import pandas as pd
import warnings
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y']

# Methods
def _update_AA_formulae(model, letter_to_met):

    # Amino acid letter to formula mapping
    amino_acid_formulas = {
        "A": "C3H5NO",
        "R": "C6H12N4O",
        "N": "C4H6N2O2",
        "D": "C4H5NO3",
        "C": "C3H5NOS",
        "Q": "C5H8N2O2",
        "E": "C5H7NO3",
        "G": "C2H3NO",
        "H": "C6H7N3O",
        "I": "C6H11NO",
        "L": "C6H11NO",
        "K": "C6H12N2O",
        "M": "C5H9NOS",
        "F": "C9H9NO",
        "P": "C5H7NO",
        "S": "C3H5NO2",
        "T": "C4H7NO2",
        "W": "C11H10N2O",
        "Y": "C9H9NO2",
        "V": "C5H9NO",
    }

    # Update each metabolite's formula in the model
    for letter, metab_id in letter_to_met.items():
        metab = model.metabolites.get_by_id(metab_id)
        metab.formula = amino_acid_formulas[letter]

    return(model)

def _get_protein_sequence(path_to_genbank):
    # Get the prot id and sequence of each protein from genbank file
    from Bio import SeqIO
    genome_record = SeqIO.parse(path_to_genbank, 'genbank')
    seq_dict = {}
    for record in genome_record:
        for i in record.features:
            if i.type == 'CDS' and 'protein_id' in i.qualifiers:
                seq_dict[i.qualifiers['protein_id'][0]] = i.qualifiers['translation'][0]

    return seq_dict

def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        return cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        return cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')

def _import_proteomic(path_to_proteomic,seq_dict):
    import pandas as pd
    import warnings
    proteomics =pd.read_csv(path_to_proteomic,header=None)
    #1- Verify number of columns
    if len(proteomics.columns) > 2:
        raise Exception("Your file format is not appropriate, more than 2 columns")
    #2- Verify presence of header
    if type(proteomics.loc[0, 0]) == str and type(proteomics.loc[0, 1]) == str:
        proteomics = proteomics.iloc[1:]
    #3- Remove null data
    if proteomics.isnull().values.any():
        proteomics = proteomics.dropna()
    #4- Verify column order (identifiers first, abundance second)
    try:
        try:
            abundances = [float(i) for i in proteomics.iloc[0:, 1]]
            identifiers = [str(i) for i in proteomics.iloc[0:, 0]]
        except:
            abundances = [float(i) for i in proteomics.iloc[0:, 0]]
            identifiers = [str(i) for i in proteomics.iloc[0:, 1]]
    except Exception:
        raise Exception('The abundances cannot be converted to float.')

    conform_df = pd.DataFrame({'identifiers':identifiers,'abundances':abundances},columns=['identifiers','abundances'])
    #5- Verify redundancy in protein identifiers
    if len(set(conform_df['identifiers'])) == len(conform_df['identifiers']):
        pass
    else:
        warnings.warn('Redundancy in dataset identifiers')
    #6- Make sure that protein id are used
    if len(list(set(conform_df['identifiers']).intersection(set(seq_dict.keys())))) > 0:
        pass
    else:
        raise Exception('Identifiers not protein_id')

    keys = [k for k in conform_df.identifiers]
    values = [v for v in conform_df.abundances]
    return dict(zip(keys, values))

def _get_aa_composition(seq_dict):
    # For each protein find the amino acid composition
    # Outputs a dictionnary of dictionnaries where:
    # Keys = locus_tag
    # Values = A dictionary for each amino acid
    # This dictionary contains:
    # Keys = amino acid by letter code
    # Values = the occurence of that amino acid
    list_of_dict = []
    for k,v in seq_dict.items():
        list_of_occurences = []
        # Get the occurence for each letter
        for letter in AMINO_ACIDS:
            protein_sequence = v
            occurence_of_letter = protein_sequence.count(letter)
            list_of_occurences.append(occurence_of_letter)
        # Generate dictionary of occurences for a given gene
        dict_of_occurences = dict(zip(AMINO_ACIDS, list_of_occurences))
        # Generate dict for each gene
        dict_per_prot = {k: dict_of_occurences}
        # Store the amount of each amino acid per gene in a list
        list_of_dict.append(dict_per_prot)

    return list_of_dict

def _normalize_aa_composition(list_of_dict, path_to_proteomic,seq_dict):
    # Normalize the value of each amino acid per protein following proteomic data
    normalized_dict = {'A': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0., 'G': 0., 'H': 0., 'I': 0.,
                       'K': 0., 'L': 0., 'M': 0., 'N': 0., 'P': 0., 'Q': 0., 'R': 0., 'S': 0., 'T': 0., 'V': 0.,
                       'W': 0., 'Y': 0.}
    # Import proteomic data into dictionnary
    proteomics = _import_proteomic(path_to_proteomic,seq_dict)
    for d in list_of_dict:
        # Get the coefficient from proteomics
        coeff = proteomics.get(str(list(d.keys())[0]))
        # If no protein abundance coefficient is 0.
        try:
            coeff_number = float(coeff)
        except:
            coeff_number = 0.

        # Multiply each amino acid by the coefficient
        amino_acids = list(d.values())
        for letter in AMINO_ACIDS:
            value = float(amino_acids[0].get(letter))
            # Update the normalized value
            normalized_value = value * coeff_number
            new_value = normalized_dict.get(letter) + normalized_value
            normalized_dict[letter] = new_value

    return normalized_dict

def _get_norm_sum(normalized_dict):
    # 1- Sum normalized ratios
    norm_sum = 0.
    for letter in AMINO_ACIDS:
        value = normalized_dict.get(letter)
        norm_sum = value + norm_sum

    return norm_sum

def _get_ratio(normalized_dict, norm_sum, PROTEIN_WEIGHT_FRACTION, CELL_WEIGHT):
    # 2- Divide letter to norm_sum to get ratio of each amino acid in the cell
    # based on proteomic data
    ratio_dict = {'A': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0., 'G': 0., 'H': 0., 'I': 0.,
                  'K': 0., 'L': 0., 'M': 0., 'N': 0., 'P': 0., 'Q': 0., 'R': 0., 'S': 0., 'T': 0., 'V': 0.,
                  'W': 0., 'Y': 0.}

    # Constant for the amount of protein in the cell
    PROTEIN_WEIGHT = CELL_WEIGHT * PROTEIN_WEIGHT_FRACTION
    for letter in AMINO_ACIDS:
        value = normalized_dict.get(letter)
        ratio = value / norm_sum
        # Convert ratios to grams
        converted_ratio = ratio * PROTEIN_WEIGHT
        ratio_dict[letter] = converted_ratio

    return ratio_dict

def _convert_to_coefficient(ratio_dict, model, letter_to_met, CELL_WEIGHT):
    WATER_WEIGHT = 18.01528

    # 3- Convert gram ratios to mmol/g Dry weight
    '''
    To verify that the normalized to grams to get to the total amount of protein
    (here everything is converted to grams instead of femto grams)
    '''

    metabolites, coefficients = [],[]
    # Get number of moles from number of grams
    for letter in AMINO_ACIDS:
        metab_id = letter_to_met.get(letter)
        metab = model.metabolites.get_by_id(metab_id)
        mol_weight = metab.formula_weight - WATER_WEIGHT
        grams = ratio_dict.get(letter)
        mmols_per_cell = (grams / mol_weight) * 1000
        mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
        coefficients.append(mmols_per_gDW)
        metabolites.append(metab)

    Protein_biomass_coefficients = dict(zip(metabolites,[-i for i in coefficients]))
    return Protein_biomass_coefficients


def generate_coefficients(path_to_genbank, 
                          path_to_model, 
                          path_to_proteomic, 
                          PROTEIN_WEIGHT_FRACTION=0.55, 
                          letter_to_met = None,
                          h2o = None):
    """

    Generates a dictionary of metabolite:coefficients for the 20 amino acids contained in proteins from the organism's
    GenBank annotated file, total Protein weight percentage and proteomic data.

    :param path_to_genbank: a path to the GenBank annotation file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param path_to_proteomic: a two column pandas dataframe (protein_id, abundance)

    :param PROTEIN_RATIO: the ratio of DNA in the entire cell

    :param letter_to_met: a dictionary with keys for each amino acid and 
        values the corresponding metabolite identifiers (string) in the model. 
        Default = None, in case which BIGG nomenclature is assumed.

    :param h2o: a string indicating the metabolite identifier for water (h2o) in the model.
        Default = None, in case which BIGG nomenclature is assumed.

    :return: a dictionary of metabolites and coefficients
    """
    model = _import_model(path_to_model)

    # Dryweight of an e. coli cell in femtogram
    # Cancels out in these equations eventually, but makes variables interpretable
    CELL_WEIGHT = 280

    if PROTEIN_WEIGHT_FRACTION > 1.:
        raise Exception('WEIGHT FRACTION should be a number between 0 and 1')
    
    if letter_to_met is None:
        letter_to_met = {'A': 'ala__L_c', 'C': 'cys__L_c',
                         'D': 'asp__L_c', 'E': 'glu__L_c',
                         'F': 'phe__L_c', 'G': 'gly_c',
                         'H': 'his__L_c', 'I': 'ile__L_c',
                         'K': 'lys__L_c', 'L': 'leu__L_c',
                         'M': 'met__L_c', 'N': 'asn__L_c',
                         'P': 'pro__L_c', 'Q': 'gln__L_c',
                         'R': 'arg__L_c', 'S': 'ser__L_c',
                         'T': 'thr__L_c', 'V': 'val__L_c',
                         'W': 'trp__L_c', 'Y': 'tyr__L_c'}
    if not all(value in [met.id for met in model.metabolites] for value in letter_to_met.values()):
        raise ValueError(f"One or more metabolite IDs from `letter_to_met` are not found in the model.")
    
    if any(model.metabolites.get_by_id(letter_to_met[letter]).formula_weight is None for letter in AMINO_ACIDS):
        warnings.warn('The molecular weight of at least one AA metabolite cannot be determined, potentially due to unspecified R-chains.' \
        'The molecular weight will be determined automatically.')
        model = _update_AA_formulae(model,letter_to_met)
        
    if h2o is None:
        h2o = 'h2o_c'
    if not h2o in [met.id for met in model.metabolites]:
        raise ValueError(f"The metabolite ID '{h2o}' is not found in the model.")

    # Operations
    # 1- Parse the genome, extract protein sequence, count and store amino acid composition of each protein
    if PROTEIN_WEIGHT_FRACTION > 1.:
        raise ValueError('Must enter ratio, value between 0. and 1.')
    seq_dict = _get_protein_sequence(path_to_genbank)
    list_of_dict = _get_aa_composition(seq_dict)
    normalized_dict = _normalize_aa_composition(list_of_dict,path_to_proteomic,seq_dict)

    # 2- Get coefficients from experimental proteomics data
    # Proteomics data should come in a 2 columns standard format protein_id:abundance
    norm_sum = _get_norm_sum(normalized_dict)
    ratio_dict = _get_ratio(normalized_dict, norm_sum, PROTEIN_WEIGHT_FRACTION, CELL_WEIGHT)
    biomass_coefficients = _convert_to_coefficient(ratio_dict,model,letter_to_met,CELL_WEIGHT)
    h2o_coeff = sum(biomass_coefficients.values())
    h2o_dict = {model.metabolites.get_by_id(h2o): -h2o_coeff}
    biomass_coefficients.update(h2o_dict)
    return biomass_coefficients

'''
The option to update the coefficients of the metabolites in the biomass objective function is left to the user
'''
def update_biomass_coefficients(dict_of_coefficients,model):
    """

    Updates the biomass coefficients given the input metabolite:coefficient dictionary.

    :param dict_of_coefficients: dictionary of metabolites and coefficients

    :param model: model to update

    :return: The biomass objective function is updated.
    """
    from BOFdat import update
    update.update_biomass(dict_of_coefficients,model)
