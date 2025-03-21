# Function to change the iso_label part of all data files for a field
#
# - written for IQFs after reordering the newforms at each level;
#   will need  a little work for other field types
#
# - assumptions:
#
# (1) we already have complete data for the field in question
#
# (2) we have rerun a curve search for that field using a newform file and the optional parameter
#     old_curves_file pointing to the existing curves file, so that no searching is done and we will
#     have exact matched of a-invariants, not just isomorphisms of curves
#
# The output will be to new files with suffix

from files import ECNF_DIR, read_all_field_data, read_curves_magma, write_data_files
from fields import nf_lookup, get_field_type_from_label
from codec import ainvs_from_string


def change_iso_label(Edata, perm, verbose=False):
    r"""
    Edata is a complete curve record, perm is a dict from old
    iso_codes to new.  The Edata is changed in place, nothing is
    returned.
    """
    from sage.databases.cremona import class_to_int
    old_iso_label = Edata['iso_label']
    pi = perm[Edata['conductor_label']]
    new_iso_label = pi[old_iso_label]
    if old_iso_label==new_iso_label:
        return
    check_labels(Edata)
    if verbose:
        print(f"Applying {pi} to\n {Edata=}:")
    for k in ['short_class_label', 'class_label', 'short_label', 'label']:
        Edata[k] = Edata[k].replace(old_iso_label, new_iso_label)
    Edata['iso_label'] = new_iso_label
    Edata['iso_nlabel'] = class_to_int(new_iso_label)
    check_labels(Edata)
    if verbose:
        print(f" -->  {Edata}")

def check_labels(Edata):
    r"""
    Given a complete curve record, check that its various label components are in place and consistent.
    """
    from sage.databases.cremona import class_to_int
    iso_label = Edata['iso_label']
    clabel = Edata['conductor_label']
    flabel = Edata['field_label']
    n = str(Edata['number'])
    assert Edata['short_class_label'] == "-".join([clabel,iso_label])
    assert Edata['class_label'] == "-".join([flabel,clabel,iso_label])
    assert Edata['short_label'] == "".join([Edata['short_class_label'], n])
    assert Edata['label'] == "".join([Edata['class_label'], n])
    assert Edata['iso_nlabel'] == class_to_int(iso_label)

def get_data_and_perm(field_label, verbose=False):
    r"""
    Given a field label, read all the old data for that field and also
    a new rawcurves file from the rawdata subdirectory of ECNF_DIR.

    Return the complete old data (list of complecte curve records) and
    perm, a dict with key=conductor_label, value=dict mapping
    old_iso_label to new_iso_label for that conductor.

    Strategy: for each curve in the rawcurves file, find a curve in
    the old data which matches (conductor_norm, conductor_label,
    ainvs) and add to the perm dictionary.

    """
    field_type = get_field_type_from_label(field_label)
    ECNF_SUBDIR = "/".join([ECNF_DIR, field_type])
    RAW_SUBDIR = "/".join([ECNF_DIR, 'rawdata'])
    rawcurves_file = RAW_SUBDIR+"/rawcurves."+field_label

    K = nf_lookup(field_label)
    data = read_all_field_data(ECNF_SUBDIR, field_label, check_cols=True)
    perm = {}
    allfound = True
    for curve in read_curves_magma(rawcurves_file):
        conductor_label = curve['conductor_label']
        if conductor_label not in perm:
            perm[conductor_label] = {}
        iso_label = curve['iso_label']
        new_label = "{}-{}-{}".format(curve['field_label'], conductor_label, iso_label)
        ainvs = curve['ainvs']
        found = False
        for old_label, Edata in data.items():
            if Edata['conductor_norm'] != curve['conductor_norm']:
                continue
            if Edata['conductor_label'] != curve['conductor_label']:
                continue
            if ainvs != ainvs_from_string(K, Edata['ainvs']):
                continue
            found = True
            if verbose:
                print(f"Curve with new label {new_label} matches old {old_label}")
            perm[conductor_label][Edata['iso_label']] = iso_label
            break
        if not found:
            print(f"Curve not found for new label {new_label}, {ainvs = }")
            allfound = False
    if allfound:
        print("all curves found")
    return data, perm

def permute_iso_labels(field_label, verbose=False):
    r"""
    For one field, read old data and new rawcurves, find the
    permutation for each conductor and apply the permutation to each
    record, returning the revised list of records.
    """
    print(f"Working on field {field_label}")
    data, perm = get_data_and_perm(field_label, verbose)
    print(" - found permutations, applying them...")
    n=0
    for old_label, Edata in data.items():
        change_iso_label(Edata, perm)
        n+=1
    print(f"done: changed {n} curve labels")
    return data

def rewrite(field_label):
    r"""
    For one field, compute the revised list of records and output to new datafiles.
    """
    data = permute_iso_labels(field_label)
    write_data_files(data,
                     field_type=get_field_type_from_label(field_label),
                     field_label=field_label,
                     base_dir=ECNF_DIR,
                     append=False,
                     suffix='reordered')




