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

# schema for SQL table ec_nfcurves:

ec_nfcurves_schema = {
    'id': 'bigint', 'torsion_structure': 'jsonb', 'ainvs': 'text', 'cm':
    'integer', 'torsion_gens': 'jsonb', 'number': 'smallint', 'rank':
    'smallint', 'field_label': 'text', 'galois_images': 'jsonb',
    'conductor_norm': 'bigint', 'abs_disc': 'bigint', 'iso_nlabel':
    'smallint', 'rank_bounds': 'jsonb', 'conductor_ideal': 'text',
    'base_change': 'jsonb', 'local_data': 'jsonb', 'analytic_rank':
    'smallint', 'minD': 'text', 'label': 'text', 'jinv': 'text',
    'conductor_label': 'text', 'class_deg': 'integer', 'reg': 'numeric',
    'class_size': 'smallint', 'isogeny_degrees': 'jsonb', 'class_label':
    'text', 'iso_label': 'text', 'degree': 'smallint', 'non_min_p':
    'jsonb', 'q_curve': 'boolean', 'short_label': 'text',
    'short_class_label': 'text', 'isogeny_matrix': 'jsonb',
    'torsion_order': 'smallint', 'non-surjective_primes': 'jsonb',
    'equation': 'text', 'gens': 'jsonb', 'ngens': 'smallint',
    'signature': 'jsonb', 'trace_hash': 'bigint', 'heights': 'numeric[]',
    'isodeg': 'integer[]', 'omega': 'numeric',
    'potential_good_reduction': 'boolean', 'tamagawa_product': 'integer',
    'bad_primes': 'jsonb', 'Lvalue': 'numeric', 'sha': 'integer',
    'torsion_primes': 'integer[]', 'n_bad_primes': 'integer',
    'reducible_primes': 'integer[]', 'semistable': 'boolean', 'disc':
    'text', 'normdisc': 'bigint',
}
