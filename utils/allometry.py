import numpy as np

def biomass_coeffs_abies(x):
    if x < 0.35:
        return (-2.3123, 2.3482)
    else:
        return (-3.1774, 2.6426)

def biomass_coeffs_cupressaceae(x):
    if x < 0.3:
        return (-1.9615, 2.1063)
    if 0.3 <= x < 0.39:
        return (-2.7765, 2.4195)
    else:
        return (-2.6327, 2.4757)

def biomass_coeffs_larix(_):
    return (-2.3012, 2.3853)

def biomass_coeffs_picea(x):
    if x < 0.35:
        return -3.0300, 2.5567
    else:
        return -2.1364, -2.3233

def biomass_coeffs_pinus(x):
    if x < 0.45:
        return -2.6177, 2.4638
    else:
        return -3.0506, 2.6465

def biomass_coeffs_pseudotsuga(_):
    return -2.4623, 2.4852

def biomass_coeffs_tsuga(x):
    if x < 0.4:
        return -2.3480, 2.3876
    else:
        return -2.9208, 2.5697

def biomass_coeffs_aceraceae(x):
    if x < 0.5:
        return -2.0470, 2.3852
    else:
        return -1.8011, 2.3852

def biomass_coeffs_betulaceae(x):
    if x < 0.4:
        return -2.5932, 2.5349
    elif 0.40 <= x <=0.49:
        return -2.2271, 2.4513
    elif 0.5 <= x <= 0.59:
        return -1.8096, 2.3480
    else:
        return -2.2652, 2.5349

def biomass_coeffs_cornaceae(_):
    return -2.2118, 2.4133

def biomass_coeffs_fabaceae_carya(_):
    return -2.5095, 2.6175

def biomass_coeffs_fabaceae_other(_):
    return -2.5095, 2.5437

def biomass_coeffs_fabaceae_deciduous(_):
    return -2.0705, 2.4410

def biomass_coeffs_fabaceae_evergreen(_):
    return -2.2198, 2.4410

def biomass_coeffs_hamamelidaceae(_):
    return -2.6390, 2.5466

def biomass_coeffs_hippocastanaceae(_):
    return -2.4108, 2.4177

def biomass_coeffs_magnoliaceae(_):
    return -2.5497, 25011

def biomass_coeffs_oleaceae(x):
    if x < 0.55:
        return -2.0314, 2.3524
    else:
        return -1.8384, 2.3524

def biomass_coeffs_salicaceae(x):
    if x < 0.35:
        return -2.6863, 2.4561
    else:
        return -2.4441, 2.4561

def cal_biomass(b1, b2, d):
    ln_biomass = b1 + b2*np.log(d)
    return np.exp(ln_biomass)

biomass_coeffs = {
    'Abies': biomass_coeffs_abies,
    'Cupressaceae': biomass_coeffs_cupressaceae,
    'Larix': biomass_coeffs_larix,
    'Picea': biomass_coeffs_picea,
    'Pinus': biomass_coeffs_pinus,
    'Pseudotsuga': biomass_coeffs_pseudotsuga,
    'Tsuga': biomass_coeffs_tsuga,
    'Aceraceae': biomass_coeffs_aceraceae,
    'Betulaceae': biomass_coeffs_betulaceae,
    'Cornaceae': biomass_coeffs_cornaceae,
    'Ericaceae': biomass_coeffs_cornaceae,
    'Lauraceae': biomass_coeffs_cornaceae,
    'Platanaceae': biomass_coeffs_cornaceae,
    'Rosaceae': biomass_coeffs_cornaceae,
    'Ulmaceae': biomass_coeffs_cornaceae,
    'Fabaceae': biomass_coeffs_fabaceae_carya,
    'Juglandaceae': biomass_coeffs_fabaceae_carya,
    'Carya': biomass_coeffs_fabaceae_carya,
    'Fabaceae_deciduous': biomass_coeffs_fabaceae_deciduous,
    'Fabaceae_evergeen': biomass_coeffs_fabaceae_evergreen,
    'Hamamelidaceae': biomass_coeffs_hamamelidaceae,
    'Hippocastanaceae': biomass_coeffs_hippocastanaceae,
    'Tiliaceae': biomass_coeffs_hippocastanaceae,
    'Magnoliaceae': biomass_coeffs_magnoliaceae,
    'Oleaceae': biomass_coeffs_oleaceae,
    'Salicaceae': biomass_coeffs_salicaceae,
}