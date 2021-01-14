# Column names for the different raw data files

label_cols =  ['field_label', 'conductor_label', 'iso_label', 'number']

column_names = { 'curves': label_cols + ['conductor_ideal', 'conductor_norm',
                                         'ainvs', 'jinv',
                                         'D', 'Dnorm',
                                         'equation',
                                         'cm', 'base_change', 'q_curve'],
                 'isoclass': label_cols + ['isogeny_matrix', 'trace_hash'],
                 'local_data': label_cols + ['local_data', 'non_min_p', 'minD'],
                 'mwdata': label_cols + ['rank', 'rank_bounds', 'analytic_rank', 'ngens', 'gens',
                                         'heights', 'reg', 'torsion_order', 'torsion_structure',
                                         'torsion_gens', 'omega', 'Lvalue', 'sha'],
                 'galrep': ['label', 'galois_images']
}
