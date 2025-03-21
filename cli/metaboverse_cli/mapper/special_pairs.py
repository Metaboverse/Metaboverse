"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) Jordan A. Berg, The University of Utah
"""

# Define redox pairs with clear separation between oxidized and reduced forms
REDOX_PAIRS = {
    # NAD+/NADH pair
    'NAD+/NADH': {
        'oxidized': ['NAD+', 'NAD', 'nicotinamide adenine dinucleotide', 'NAD (oxidized)', 'β-nicotinamide adenine dinucleotide'],
        'reduced': ['NADH', 'nicotinamide adenine dinucleotide reduced', 'NADH (reduced)', 'β-nicotinamide adenine dinucleotide reduced']
    },
    
    # NADP+/NADPH pair
    'NADP+/NADPH': {
        'oxidized': ['NADP+', 'NADP', 'nicotinamide adenine dinucleotide phosphate', 'NADP (oxidized)', 'β-nicotinamide adenine dinucleotide phosphate'],
        'reduced': ['NADPH', 'nicotinamide adenine dinucleotide phosphate reduced', 'NADPH (reduced)', 'β-nicotinamide adenine dinucleotide phosphate reduced']
    },
    
    # FAD/FADH2 pair
    'FAD/FADH2': {
        'oxidized': ['FAD', 'flavin adenine dinucleotide', 'flavin adenine dinucleotide oxidized', 'FAD (oxidized)'],
        'reduced': ['FADH2', 'flavin adenine dinucleotide reduced', 'FADH2 (reduced)']
    },
    
    # CoQ/CoQH2 pair
    'CoQ/CoQH2': {
        'oxidized': ['CoQ', 'ubiquinone', 'coenzyme Q', 'coenzyme Q10', 'CoQ10', 'ubiquinone-10'],
        'reduced': ['CoQH2', 'ubiquinol', 'coenzyme QH2', 'reduced coenzyme Q', 'reduced coenzyme Q10', 'CoQ10H2', 'ubiquinol-10']
    },
    
    # GSH/GSSG pair
    'GSH/GSSG': {
        'oxidized': ['GSSG', 'glutathione disulfide', 'oxidized glutathione'],
        'reduced': ['GSH', 'glutathione', 'reduced glutathione', 'L-glutathione']
    },
    
    # Trx-S2/Trx-(SH)2 pair
    'Trx-S2/Trx-(SH)2': {
        'oxidized': ['Trx-S2', 'thioredoxin oxidized', 'oxidized thioredoxin'],
        'reduced': ['Trx-(SH)2', 'thioredoxin reduced', 'reduced thioredoxin']
    },
    
    # Cytochrome C pair
    'Cytochrome C (Fe3+)/Cytochrome C (Fe2+)': {
        'oxidized': ['cytochrome c (Fe3+)', 'cytochrome c oxidized', 'oxidized cytochrome c'],
        'reduced': ['cytochrome c (Fe2+)', 'cytochrome c reduced', 'reduced cytochrome c']
    },
    
    # Ferredoxin pair
    'Ferredoxin (oxidized)/Ferredoxin (reduced)': {
        'oxidized': ['ferredoxin oxidized', 'oxidized ferredoxin'],
        'reduced': ['ferredoxin reduced', 'reduced ferredoxin']
    },
    
    # Ascorbic acid pair
    'Ascorbic acid/Dehydroascorbic acid': {
        'oxidized': ['dehydroascorbic acid', 'DHA', 'oxidized vitamin C'],
        'reduced': ['ascorbic acid', 'vitamin C', 'L-ascorbic acid', 'ascorbate']
    },
    
    # Lipoic acid pair
    'Lipoic acid/Dihydrolipoic acid': {
        'oxidized': ['lipoic acid', 'α-lipoic acid', 'thioctic acid', 'oxidized lipoic acid'],
        'reduced': ['dihydrolipoic acid', 'DHLA', 'reduced lipoic acid']
    },
    
    # FMN/FMNH2 pair
    'FMN/FMNH2': {
        'oxidized': ['FMN', 'flavin mononucleotide', 'riboflavin 5\'-phosphate', 'oxidized FMN'],
        'reduced': ['FMNH2', 'reduced flavin mononucleotide', 'reduced riboflavin 5\'-phosphate', 'reduced FMN']
    }
}

def get_redox_state(metabolite_name, synonyms):
    """Determine if a metabolite is part of a redox pair and its state
    
    Args:
        metabolite_name (str): The name of the metabolite to check
        synonyms (list): List of synonyms for the metabolite
        
    Returns:
        tuple: (redox_pair_name, redox_state) where redox_state is 'oxidized' or 'reduced'
               Returns (None, None) if not a redox pair
    """
    metabolite_lower = metabolite_name.lower()
    
    for pair_name, forms in REDOX_PAIRS.items():
        # Check if the metabolite name matches any of the forms
        if metabolite_lower in [s.lower() for s in forms['oxidized']]:
            return pair_name, 'oxidized'
        if metabolite_lower in [s.lower() for s in forms['reduced']]:
            return pair_name, 'reduced'
            
        # Check synonyms if no direct match
        for syn in synonyms:
            syn_lower = syn.lower()
            if syn_lower in [s.lower() for s in forms['oxidized']]:
                return pair_name, 'oxidized'
            if syn_lower in [s.lower() for s in forms['reduced']]:
                return pair_name, 'reduced'
    
    return None, None

def filter_redox_synonyms(synonyms, redox_pair, redox_state):
    """Filter synonyms to only include those matching the correct redox state
    
    Args:
        synonyms (list): List of synonyms to filter
        redox_pair (str): Name of the redox pair (e.g. 'NAD+/NADH')
        redox_state (str): Either 'oxidized' or 'reduced'
        
    Returns:
        list: Filtered list of synonyms that match the correct redox state
    """
    if not redox_pair or not redox_state:
        return synonyms
        
    allowed_synonyms = set(REDOX_PAIRS[redox_pair][redox_state])
    return [s for s in synonyms if s.lower() in [a.lower() for a in allowed_synonyms]]