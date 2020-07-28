def to_mass_frac(mol, vmr, mu):
    from taurex.util.util import calculate_weight
    return calculate_weight(mol)*vmr/mu