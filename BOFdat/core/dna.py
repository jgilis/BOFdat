"""
DNA
===

This module generates BOFsc for the 4 bases of DNA (dATP, dTTP, dCTP and dGTP)

"""
import warnings
from collections import Counter

BASES = ['A', 'T', 'C', 'G']

# Methods
def _import_genome(fasta):
    from Bio import SeqIO
    try:
        #Import as a single handle genome    
        genome = list(SeqIO.parse(fasta,'fasta'))
        if len(genome) > 1:
            warnings.warn('%s handles in the genome file.This may indicate that your genome is not completely assembled. \nBOFdat will parse the contigs but the stoichiometric coefficients may not be accurate.'%(len(genome),))
    except:
        raise ImportError('The file provided cannot be imported.')
        
    return genome

def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        return cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        return cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')

def _get_number_of_bases(genome, ignoreLower=False):
    total_counts = Counter()
    for record in genome:
        total_counts.update(str(record.seq))
    
    if ignoreLower:
        for base in ('a', 't', 'c', 'g'):
            total_counts.pop(base, None)
    else:
        for base_upper, base_lower in [('A', 'a'), ('T', 't'), ('C', 'c'), ('G', 'g')]:
            total_counts[base_upper] += total_counts.pop(base_lower, 0)

    return {base: total_counts.get(base, 0) for base in ['A', 'T', 'C', 'G']}

def _get_ratio(base_genome):
    # Get the ratios for each letter in the genome
    # TODO running the original _get_ratio leads to ratios that sum to 0.5, not 1
    # While this may be intentional, it is very strange.
    # Below, I adapted the code to what I think it should be.
    # This will be equivalent to running the original code and 
    # multiplying the ratios, or the downstream coefficients, by 2
    total = 2*sum(base_genome.values()) # double-stranded
    at_sum = base_genome['A'] + base_genome['T']
    gc_sum = base_genome['G'] + base_genome['C']
    return {k: at_sum / total if k in 'AT' else gc_sum / total for k in 'ATGC'}

def _convert_to_coefficient(model, ratio_genome, DNA_RATIO, base_to_met):
    
    DIPHOSPHATE_WEIGHT = 174.951262

    # Dryweight of an e. coli cell in femtogram
    # Cancels out in these equations, but makes variables interpretable
    CELL_WEIGHT = 280

    # Transform the ratios into mmol/gDW
    DNA_WEIGHT = CELL_WEIGHT * DNA_RATIO

    if base_to_met is None:
        base_to_met = {'A': 'datp_c', 'T': 'dttp_c',
                       'C': 'dctp_c', 'G': 'dgtp_c'}

    coefficients,metabolites = [],[]

    # Calculate the biomass coefficient for each metabolite
    for letter in BASES:
        ratio = ratio_genome.get(letter)
        total_weight = ratio * DNA_WEIGHT
        metab_id = base_to_met.get(letter)
        metab = model.metabolites.get_by_id(metab_id)
        mol_weight = metab.formula_weight - DIPHOSPHATE_WEIGHT
        mmols_per_cell = (total_weight / mol_weight) * 1000
        mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
        coefficients.append(mmols_per_gDW)
        metabolites.append(base_to_met.get(letter))

    DNA_coefficients = dict(zip(metabolites,[-i for i in coefficients]))
    return DNA_coefficients

def generate_coefficients(path_to_fasta,
                          path_to_model ,
                          DNA_WEIGHT_FRACTION=0.031,
                          base_to_met = None,
                          ppi = None):
    """
    Generates a dictionary of metabolite:coefficients for the 4 DNA bases from the organism's
    DNA fasta file and the weight percentage of DNA in the cell.

    :param path_to_fasta: a path to the DNA fasta file of the organism, format should be compatible 
    with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param DNA_RATIO: the ratio of DNA in the entire cell

    :param base_to_met: a dictionary with keys for each nucleotide (ACGT) and 
    values the corresponding metabolite identifiers (string) in the model. 
    Default = None, in case which BIGG nomenclature is assumed.

    :param ppi: a string indicating the metabolite identifier for pyrophosphate
    (ppi) in the model.
    Default = None, in case which BIGG nomenclature is assumed.

    :return: a dictionary of metabolites and coefficients
    """
    if DNA_WEIGHT_FRACTION > 1.:
        raise Exception('WEIGHT FRACTION should be a number between 0 and 1')
    #Operations
    genome = _import_genome(path_to_fasta)
    base_in_genome = _get_number_of_bases(genome)
    ratio_in_genome = _get_ratio(base_in_genome)
    model = _import_model(path_to_model)
    biomass_coefficients = _convert_to_coefficient(model,
                                                   ratio_in_genome,
                                                   DNA_WEIGHT_FRACTION, 
                                                   base_to_met)
    # Add pyrophosphate synthesis as the sum of the coefficients
    ppi_coeff = sum(biomass_coefficients.values())
    if ppi is None:
        ppi_dict = {model.metabolites.get_by_id('ppi_c'):-ppi_coeff}
    else:
        ppi_dict = {model.metabolites.get_by_id(ppi):-ppi_coeff}
    biomass_coefficients.update(ppi_dict)

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
