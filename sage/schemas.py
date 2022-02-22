from six import text_type
from sage.all import Set, ZZ

# Column names for the different raw data files

label_cols =  ['field_label', 'conductor_label', 'iso_label', 'number']

column_names = { 'curves': label_cols + ['conductor_ideal', 'conductor_norm',
                                         'ainvs', 'jinv',
                                         'disc', 'normdisc',
                                         'equation',
                                         'cm', 'base_change', 'q_curve'],
                 'isoclass': label_cols + ['isogeny_matrix', 'trace_hash'],
                 'local_data': label_cols + ['local_data', 'non_min_p', 'minD'],
                 'mwdata': label_cols + ['rank', 'rank_bounds', 'analytic_rank', 'ngens', 'gens',
                                         'heights', 'reg', 'torsion_order', 'torsion_structure',
                                         'torsion_gens', 'omega', 'Lvalue', 'sha'],
                 'galrep': ['label', 'galois_images']
}

all_file_types = list(column_names.keys())

# schema for SQL table ec_nfcurves:

ec_nfcurves_schema = {
    'id': 'bigint', 'torsion_structure': 'jsonb', 'ainvs': 'text', 'cm':
    'integer', 'cm_type': 'smallint', 'torsion_gens': 'jsonb', 'number': 'smallint', 'rank':
    'smallint', 'field_label': 'text', 'galois_images': 'jsonb',
    'conductor_norm': 'bigint', 'abs_disc': 'bigint', 'iso_nlabel':
    'smallint', 'rank_bounds': 'jsonb', 'conductor_ideal': 'text',
    'base_change': 'jsonb', 'local_data': 'jsonb', 'analytic_rank':
    'smallint', 'minD': 'text', 'label': 'text', 'jinv': 'text',
    'conductor_label': 'text', 'class_deg': 'integer', 'reg': 'numeric',
    'class_size': 'smallint', 'class_label':
    'text', 'iso_label': 'text', 'degree': 'smallint', 'non_min_p':
    'jsonb', 'q_curve': 'boolean', 'short_label': 'text',
    'short_class_label': 'text', 'isogeny_matrix': 'jsonb',
    'torsion_order': 'smallint', 'nonmax_primes': 'smallint[]', 'nonmax_rad': 'integer',
    'equation': 'text', 'gens': 'jsonb', 'ngens': 'smallint',
    'signature': 'jsonb', 'trace_hash': 'bigint', 'heights': 'numeric[]',
    'isodeg': 'integer[]', 'omega': 'numeric',
    'potential_good_reduction': 'boolean', 'tamagawa_product': 'integer',
    'bad_primes': 'jsonb', 'Lvalue': 'numeric', 'sha': 'integer',
    'torsion_primes': 'integer[]', 'n_bad_primes': 'integer',
    'reducible_primes': 'integer[]', 'semistable': 'boolean', 'disc':
    'text', 'normdisc': 'numeric', 'conductor_norm_factors': 'integer[]',
    'non-surjective_primes': 'jsonb',
    'root_analytic_conductor': 'double precision',
}

# Python types of the columns in ec_nfcurves:

str_type = text_type
int_type = type(int(1))
float_type = type(float(1))
list_type = type([1,2,3])
bool_type = type(True)
hash_type = type(ZZ(2**65).__int__())

keys_and_types = {'field_label':  str_type,
                  'degree': int_type,
                  'signature': list_type, # of ints
                  'abs_disc': int_type,
                  'label':  str_type,
                  'short_label':  str_type,
                  'class_label':  str_type,
                  'short_class_label':  str_type,
                  'class_deg':  int_type,
                  'class_size':  int_type,
                  'conductor_label': str_type,
                  'conductor_ideal': str_type,
                  'conductor_norm': int_type,
                  'iso_label': str_type,
                  'iso_nlabel': int_type,
                  'number': int_type,
                  'ainvs': str_type,
                  'jinv': str_type,
                  'cm': int_type,
                  'cm_type': int_type,
                  'ngens': int_type,
                  'rank': int_type,
                  'rank_bounds': list_type, # 2 ints
                  'analytic_rank': int_type,
                  'torsion_order': int_type,
                  'torsion_structure': list_type, # 0,1,2 ints
                  'gens': list_type, # of strings
                  'torsion_gens': list_type, # of strings
                  'isogeny_matrix': list_type, # of lists of ints
                  'isodeg': list_type, # of ints
                  'class_deg': int_type,
                  'non-surjective_primes': list_type, # of ints
                  'nonmax_primes': list_type, # of ints
                  'nonmax_rad': int_type,
                  'galois_images': list_type, # of strings
                  'equation': str_type,
                  'local_data': list_type, # of dicts
                  'non_min_p': list_type, # of strings
                  'minD': str_type,
                  'disc': int_type,
                  'normdisc': str_type,
                  'heights': list_type, # of floats
                  'reg': float_type, # or int(1)
                  'q_curve': bool_type,
                  'base_change': list_type, # of strings
                  'trace_hash': hash_type,
                  'conductor_norm_factors': list_type,
                  'root_analytic_conductor': float_type,
}

ec_nfcurves_extra_columns = ['omega', 'potential_good_reduction', 'semistable', 'tamagawa_product', 'bad_primes', 'Lvalue', 'sha', 'torsion_primes', 'n_bad_primes', 'reducible_primes']

extra_keys_and_types = {'omega': float_type,
                        'potential_good_reduction': bool_type,
                        'semistable': bool_type,
                        'tamagawa_product': int_type,
                        'bad_primes': list_type,
                        'Lvalue': float_type,
                        'sha': int_type,
                        'torsion_primes': list_type,
                        'n_bad_primes': int_type,
                        'reducible_primes': list_type}

keys_and_types.update(extra_keys_and_types)

extra_keys_and_postgres_types = {'omega': 'numeric',
                        'potential_good_reduction': 'boolean',
                        'semistable': 'boolean',
                        'tamagawa_product': 'integer',
                        'bad_primes': 'jsonb',
                        'Lvalue': 'numeric',
                        'sha': 'integer',
                        'torsion_primes': 'integer[]',
                        'n_bad_primes': 'integer',
                        'reducible_primes': 'integer[]'}


ec_nfcurves_columns = ec_nfcurves_all_columns = Set(keys_and_types.keys()) + Set(['id'])
