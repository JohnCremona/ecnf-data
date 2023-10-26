# Functions for reading and writing data files

import os
import sys

from codec import (curve_from_string, curve_from_strings,
                   ainvs_from_strings, ainvs_from_string,
                   convert_ideal_label, decode_points_one2many,
                   decode_int_list,
                   encode_points, file_line, parse_curves_line,
                   parse_isoclass_line, parse_local_data_line,
                   parse_mwdata_line, parse_new_mwdata_line,
                   parse_galrep_line, parse_line_label_cols)

from fields import nf_lookup
from schemas import all_file_types
from modularitycheck import isModular

HOME = os.getenv("HOME")
BIANCHI_DATA_DIR = os.path.join(HOME, "bianchi-data")
ECNF_DIR = os.path.join(HOME, "ecnf-data")
ECNF_UPLOAD_DIR = os.path.join(HOME, "ecnf-upload")

with open(os.path.join(ECNF_DIR, "field_types")) as ftypes:
    all_ftypes = [ft.replace("\n","") for ft in ftypes.readlines()]

def read_all_field_data(base_dir, field_label, check_cols=True, mwdata_format='new'):
    r"""Given a field label, read all the data in files curves.field_label,
    isoclass.field_label, local_data.field_label, mwdata.field_label,
    galrep.field_label (from directory base_dir), returning a single
    dict with keys full labels and values complete curve records.
    Here field_label does not have to actually be a field label
    (e.g. for testing).

    If we knew that the data in all 5 files was in the same order
    (taking into account the fact that isoclass.* files only have one
    line per isogeny class) this could be a generator yielding one
    curve record at a time, but we do not want to assume that, so we
    collect everything and return it all at once.

    If check_cols is True then we check that every record has keys
    exactly as in the standard list from the table ec_nfcurves
    (omitting the 'id' column).

    """
    from schemas import ec_nfcurves_all_columns
    from sage.all import is_prime, Set
    curves_filename = os.path.join(base_dir, 'curves.{}'.format(field_label))
    isoclass_filename = os.path.join(base_dir, 'isoclass.{}'.format(field_label))
    local_data_filename = os.path.join(base_dir, 'local_data.{}'.format(field_label))
    mwdata_filename = os.path.join(base_dir, 'mwdata.{}'.format(field_label))
    galrep_filename = os.path.join(base_dir, 'galrep.{}'.format(field_label))

    all_data = {}
    n = 0

    print("Reading from {}".format(curves_filename))
    with open(curves_filename) as curves:
        for L in curves:
            label, record = parse_curves_line(L)
            if label:
                n += 1
                all_data[label] = record
    print("Read {} curves from {}".format(n,curves_filename))
    n = 0

    print("Reading from {}".format(isoclass_filename))
    with open(isoclass_filename) as isoclass:
        for L in isoclass:
            label, record = parse_isoclass_line(L)
            if label:
                n += 1
                nc = record['class_size']
                all_iso_degs = record.pop('all_iso_degs')
                # reducible_primes is the same for all curves in the class
                reducible_primes = [ d for d in all_iso_degs[1] if is_prime(d)]
                label = label[:-1] # delete final '1'
                for ic in range(1,nc+1):
                    record['reducible_primes'] = reducible_primes
                    record['isodeg'] = all_iso_degs[ic]
                    all_data[label+str(ic)].update(record)
    print("Read {} classes from {}".format(n,isoclass_filename))
    n = 0

    print("Reading from {}".format(local_data_filename))
    with open(local_data_filename) as ld:
        for L in ld:
            label, record = parse_local_data_line(L)
            if label:
                n += 1
                all_data[label].update(record)
    print("Read {} local_data records from {}".format(n,local_data_filename))
    n = 0

    print("Reading from {}".format(mwdata_filename))
    with open(mwdata_filename) as mwdata:
        for L in mwdata:
            if mwdata_format=='old':
                label, record = parse_mwdata_line(L)
            else:
                label, record = parse_new_mwdata_line(L)
            if label:
                n += 1
                all_data[label].update(record)
    print("Read {} mwdata records from {}".format(n,mwdata_filename))
    n = 0

    print("Reading from {}".format(galrep_filename))
    with open(galrep_filename) as galrep:
        for L in galrep:
            label, record = parse_galrep_line(L)
            if label:
                n += 1
                all_data[label].update(record)
    print("Read {} galrep records from {}".format(n,galrep_filename))

    if check_cols:
        for label in all_data:
            cols = Set(all_data[label]) + Set(['id'])
            if cols != ec_nfcurves_all_columns:
                print("Wrong key set for {}".format(label))
                diff = cols - ec_nfcurves_all_columns
                if diff:
                    print("data has extra keys (not in ec_nfcurves_all_columns) {}".format(diff))
                diff = ec_nfcurves_all_columns - cols
                if diff:
                    print("data is missing keys (in ec_nfcurves_all_columns) {}".format(diff))
    return all_data

def read_curve_file(infile):
    r""" Read an old-format curves file, each line containing 13 data
    fields as defined the in the ecnf-format.txt file. Output yields
    dicts, one per curve, in the same order as input.
    """
    index = 0
    with open(infile) as file:
        for L in file.readlines():
            if L[0]=='#': # allow for comment lines
                continue
            data = L.split()
            if len(data)!=13:
                print("line {} does not have 13 fields, skipping".format(L))
                continue
            index += 1
            curve = {'index': index,
                     'field_label': data[0],
                     'conductor_label': data[1],
                     'iso_label': data[2],
                     'c_num': data[3],
                     'conductor_ideal': data[4],
                     'conductor_norm': data[5],
                     'ainvs': data[6:11],
                     'cm_flag': data[11],
                     'q_curve_flag': data[12]
            }
            yield curve
    return

def read_isoclass_file(infile):
    r"""
    Read an isoclass file, each line containing 5 data fields as defined
    the in the ecnf-format.txt file. Output is a list of dicts, one
    per curve, in the same order as input.
    """
    index = 0
    with open(infile) as file:
        for L in file.readlines():
            if L[0]=='#': # allow for comment lines
                continue
            data = L.split()
            if len(data) != 5:
                print("line {} does not have 5 fields, skipping".format(L))
                continue
            index += 1
            curve = {'index': index,
                     'field_label': data[0],
                     'conductor_label': data[1],
                     'iso_label': data[2],
                     'c_num': data[3],
                     'isomat': data[4]
            }
            yield curve
    return

def read_curves(infile, only_one=False, ncurves=0):
    r"""
    Iterator to loop through lines of a curves.* file each
    containing 13 data fields as defined the in the ecnf-format.txt file,
    yielding its curves as EllipticCurve objects and other data:
       (field_label,conductor_label,conductor_ideal,iso_label,c_num,E)

    If only_one is True, skips curves whose 4th data field is
    *not* 1, hence only yielding one curve per isogeny class.
    """
    count=0
    with open(infile) as file:
        for L in file.readlines():
            data = L.split()
            if len(data)!=14:
                print("line {} does not have 13 fields, skipping".format(L))
                continue
            if only_one and data[3]!='1':
                continue
            count +=1
            if ncurves and count>ncurves:
                return
            field_label = data[0]
            K = nf_lookup(field_label).change_names('w')
            conductor_label = data[1]
            iso_label = data[2]
            c_num = data[3]
            conductor_ideal = data[4]
            E = curve_from_string(K, data[6])
            yield (field_label,conductor_label,conductor_ideal,iso_label,c_num,E)

def read_curves_new(infile, only_one=False, ncurves=0):
    r""" Iterator to loop through lines of a curves.* file each containing
    14 data fields as defined in the ecnf-format.txt file (see
    schemas.py), yielding its curves as EllipticCurve objects.

    If only_one is True, skips curves whose 4th data field is
    *not* 1, hence only yielding one curve per isogeny class.

    Output: yields records with keys

    field_label, conductor_norm, conductor_label, conductor_ideal, iso_label, ainvs
    """
    from schemas import column_names
    cols = column_names['curves']
    ncols = len(cols)
    count=0
    with open(infile) as file:
        for L in file.readlines():
            data = L.split()
            if len(data)!=ncols:
                print("line {} does not have 12 fields, skipping".format(L))
                continue
            if only_one and data[3]!='1':
                continue
            count +=1
            if ncurves and count>ncurves:
                return
            record = {}
            record['field_label'] = field_label = data[0]
            K = nf_lookup(field_label)
            K.__gens_dict().update({'w':K.gen()})
            record['conductor_label'] = data[1]
            record['iso_label'] = data[2]
            record['number'] = data[3]
            record['conductor_ideal'] = data[4]
            record['conductor_norm'] = data[5]
            record['ainvs'] = ainvs_from_string(K, data[6])
            yield record

def read_curves_magma(infile, min_norm=1, max_norm=None):
    r"""Iterator to loop through lines of a file containing output from a
    Magma search.  (Nothing in this function really relates to Magma.)
    For each curve there are 3 or 4 lines in the file, with prefixes

    Field (only needed once at the top if all curves are over the same field)
    Conductor
    Isogeny_class
    Curve

    containing, respectively:

    field_label, e.g. 6.6.980125.1
    conductor_ideal (string containing conductor ideal, e.g. [379,379,-w^5-w^4+6*w^3+5*w^2-7*w-4]
    "-".join(conductor_label, iso_label), e.g. 379.1-a
    ainvs, e.g. [w^2,-w^5+w^3-w^2+w+1,w^2+1,-2*w^5+24*w^4-13*w^3-54*w^2+32*w+6,51*w^4-91*w^3-29*w^2+45*w+7]

    where field elements are expressed in terms of a field generator
    'w'.  The field label must be one for which nf_lookup() can return
    the correct field.

    Output: yields records with keys

    field_label, conductor_norm, conductor_label, conductor_ideal, iso_label, ainvs

    (omitting any whose conductor_norm is not within the bounds given, if any)

    NB If the search failed to find a curve with given conductor,
    instead of a line starting "Curve " there will be a line "No curve
    found".

    """
    from sage.all import EllipticCurve, ZZ
    record = {}
    with open(infile) as file:
        for L in file.readlines():
            data = L.strip().split(maxsplit=1)
            if len(data) == 0:
                continue
            assert len(data) == 2, "Line '{}' does not have two fields".format(L)
            if data[0] == 'Field':
                field_label = data[1]
                # reset the record for a new curve
                record = {'field_label': field_label}
                K = nf_lookup(field_label)
                K.__gens_dict().update({'w':K.gen()})
            elif data[0] == 'Conductor':
                record['conductor_ideal'] = data[1]
                #cond = data[1].replace("[", "").replace("]", "").replace("(", "").replace(")", "")
                #N = K.ideal([K(a) for a in cond.split(",")])
            elif data[0] == 'Isogeny_class':
                conductor_label, iso_label = data[1].split("-")
                record['conductor_norm'] = ZZ(conductor_label.split(".")[0])
                record['conductor_label'] = conductor_label
                record['iso_label'] = iso_label
            elif data[0] == 'Curve':
                ainvs = data[1]
                ainvs = ainvs[1:-1].split(",")
                record['ainvs'] = ainvs = [K(ai) for ai in ainvs]
                E = EllipticCurve(K, ainvs)
                EN = E.conductor()
                N_norm = record['conductor_norm']
                assert EN.norm() == N_norm
                if N_norm >= min_norm and (max_norm is None or N_norm <= max_norm):
                    yield record
            elif data[0] == 'No':
                print(f"No curve for class {record['conductor_label']}{record['iso_label']}, skipping")
                continue
            else:
                print("Unrecognised line prefix {}, skipping this line".format(data[0]))
                continue

def read_classes(infile):
    r"""
    Iterator to loop through lines of a curves.* file each
    containing 13 data fields as defined the in the ecnf-format.txt file,
    yielding complete isogeny classes.  These are dicts with keys

    field_label, conductor_label, conductor_ideal, conductor_norm, iso_label, curves

    the last being a list of EllipticCurves.
    """
    count=0
    prev_class_id = ''
    this_class = {}
    with open(infile) as file:
        for L in file.readlines():
            data = L.split()
            if len(data)!=13:
                print("line {} does not have 13 fields, skipping".format(L))
                continue
            count +=1
            field_label = data[0]
            K = nf_lookup(field_label)
            conductor_label = data[1]
            iso_label = data[2]
            #c_num = data[3]
            conductor_ideal = data[4]
            conductor_norm = int(data[5])
            E = curve_from_strings(K, data[6:11])
            this_class_id = "-".join(data[:3])
            if this_class_id == prev_class_id:
                this_class['curves'].append(E)
            else:
                if this_class: # else we have just started
                    yield this_class
                this_class = {}
                this_class['field_label'] = field_label
                this_class['conductor_label'] = conductor_label
                this_class['conductor_ideal'] = conductor_ideal
                this_class['conductor_norm'] = conductor_norm
                this_class['iso_label'] = iso_label
                this_class['curves'] = [E]
                prev_class_id = this_class_id
    # final class:
    yield this_class

def read_classes_new(infile):
    r"""
    Iterator to loop through lines of a curves.* file each
    containing 13 data fields as defined the in the ecnf-format.txt file,
    yielding complete isogeny classes.  These are dicts with keys

    field_label, conductor_label, conductor_ideal, conductor_norm, iso_label, curves

    the last being a list of EllipticCurves.
    """
    count=0
    prev_class_id = ''
    this_class = {}
    with open(infile) as file:
        for L in file.readlines():
            count +=1
            label, record = parse_curves_line(L)
            field_label = record['field_label']
            K = nf_lookup(field_label)
            conductor_label = record['conductor_label']
            iso_label = record['iso_label']
            conductor_ideal = record['conductor_ideal']
            conductor_norm = record['conductor_norm']
            E = curve_from_string(K, record['ainvs'])
            this_class_id = record['class_label']
            if this_class_id == prev_class_id:
                this_class['curves'].append(E)
            else:
                if this_class: # else we have just started
                    yield this_class
                this_class = {}
                this_class['field_label'] = field_label
                this_class['conductor_label'] = conductor_label
                this_class['conductor_ideal'] = conductor_ideal
                this_class['conductor_norm'] = conductor_norm
                this_class['iso_label'] = iso_label
                this_class['curves'] = [E]
                prev_class_id = this_class_id
    # final class:
    yield this_class

def bmf_filename_from_D(absD):
    d=absD if absD%4 else absD//4
    return HOME + "/bianchi-data/nflist/nflist.{}.1-20000".format(d)

def read_newform_data(bmf_filename, verbose=False):
    r"""
    Returns a dict containing Bianchi newforms read from a file.

    INPUT:

    - ``bmf_filename`` (string) -- name of file containing Bianchi
      newform data.  Either starts with "nflist" or "newforms"; the
      latter has additional fields which we ignore here.

    - ``verbose`` (boolean, default False) -- verbosity flag.

    OUTPUT:

    dict with keys level labels (strings), values dicts with keys
    newforms labels (strings 'a', 'b', etc), values dicts with keys
    'gen_str' (string representing generator of the level as a
    principal ideal), 'sign' (sign of functional equation of the
    L-function), 'aq' (list of Atkin-Lehner eigenvalues), 'ap' (list
    of L-function coefficients of prime ideals in standard order).

    """
    bmf_file = open(bmf_filename)
    old_fmt =  "nflist" in bmf_filename
    print("file has {} format".format('old' if old_fmt else 'new'))
    newforms = {}
    for L in bmf_file.readlines():
        if verbose:
            print("raw input: %s" % L)
        if old_fmt:
            label, gen, sfe, loverp, ALs, aplist = L.split()
            level, letter = label.split("-")
        else:
            field_label, level, letter, gen, wt, bc, cm, sfe, loverp, ALs, poly, aplist = L.split()

        if not level in newforms:
            newforms[level] = {}
        nfs = newforms[level]
        nfs[letter] = nf = {}
        nf['gen_str'] = gen[1:-1]
        nf['sign'] = int(sfe)
        nf['aq'] = decode_int_list(ALs)
        nf['ap'] = decode_int_list(aplist)
        if verbose:
            print("newform data: %s" % L)

    bmf_file.close()
    return newforms

def read_missing_levels(infile):
    r"""
    Yields level labels from a file of newform/isogeny class labels, without repeats.
    """
    levels = []
    for L in infile.readlines():
        if "nflist" in infile:
            level = L.split("-")[0]
        else:
            level = L.split()[1]
        if not level in levels:
            levels += [level]
            yield level

def label_conversion_table(infile, outfile):
    with open(outfile, mode='w') as OF, open(os.path.join(BIANCHI_DATA_DIR, "ideals", infile)) as IF:
        for L in IF:
            field, ideal = L.split()
            label = convert_ideal_label(nf_lookup(field),ideal)
            OF.write(' '.join([field, ideal, label])+'\n')

def extend_mwdata(base_dir, field_label, suffix='x', minN=None, maxN=None, one_label=None, CM_only=False, max_sat_prime = None, prec=None, verbose=False):
    r"""
    Reads curves and local data files.
    Computes analytic rank and L-value using Magma, and omega (global period).
    Computes analytic Sha (rounded).
    Rewrites mwdata, inserting analytic rank value and adding columns omega, Lvalue, sha.

    The prec parameter affects the precision to which the L-value and
    global period is computed.  It is bit precision.  Magma's default
    is 6dp or 20 bits for the L-value and the running time increases
    rapidly.
    """
    from nfscripts import extend_mwdata_one
    data = read_all_field_data(base_dir, field_label, check_cols=False)
    classdata = {} # will hold isogeny-invariant values keyed by class label
    Kfactors = {} # BSD factor depending only on the field K

    from magma import get_magma

    if one_label:
        mwoutfile = base_dir+'/mwdata.'+suffix+"."+one_label
    else:
        mwoutfile = base_dir+'/mwdata.'+field_label+suffix
    with open(mwoutfile, 'w', 1) as mwdata:
        for label, Edata in data.items():
            class_label = Edata['class_label']
            if one_label and one_label!=class_label:
                continue
            if CM_only and not 'CM' in class_label:
                continue
            N = Edata['conductor_norm']
            if (minN and N<minN) or (maxN and N>maxN):
                continue
            if verbose:
                print("Processing {}".format(label))
                #print(Edata)

            if not class_label in classdata: # we'll use magma in extend_mwdata_one
                magma = get_magma()

            Edata = extend_mwdata_one(Edata, classdata, Kfactors, magma,
                                      max_sat_prime = max_sat_prime, prec=prec, verbose=verbose)
            line = file_line('mwdata', Edata)
            if verbose:
                print("New mwdata line: {}".format(line))
            mwdata.write(line + "\n")

def fix_torsion(base_dir, field_label, suffix='x', minN=None, maxN=None, one_label=None, verbose=False):
    r"""
    Reads all data files.
    Recomputes torsion data only (as this was corrupted in some cases).
    Rewrites mwdata, with tosion fields corrected.
    """
    data = read_all_field_data(base_dir, field_label, check_cols=False, mwdata_format='new')
    if one_label:
        mwoutfile = base_dir+'/mwdata.'+suffix+"."+one_label
    else:
        mwoutfile = base_dir+'/mwdata.'+field_label+suffix
    with open(mwoutfile, 'w', 1) as mwdata:
        for label, Edata in data.items():
            if one_label and one_label!=Edata['class_label']:
                continue
            N = Edata['conductor_norm']
            if (minN and N<minN) or (maxN and N>maxN):
                continue
            if verbose:
                print("Processing {}".format(label))
                #print("Edata = {}".format(Edata))
            K = nf_lookup(Edata['field_label'])
            E = curve_from_string(K,Edata['ainvs'])
            T = E.torsion_subgroup()
            nt = T.order()
            Edata['torsion_order'] = nt
            Edata['torsion_structure'] = list(T.invariants())
            tgens = [P.element() for P in T.gens()]
            Edata['torsion_gens'] = decode_points_one2many(encode_points(tgens)) # list of strings
            if verbose:
                print("Torsion order = {}".format(nt))
                print("Torsion gens = {}".format(Edata['torsion_gens']))

            line = file_line('mwdata', Edata)
            if verbose:
                print("New mwdata line: {}".format(line))
            mwdata.write(line + "\n")

def read_ai(curvefile, only_one=False, ncurves=0):
    r"""Iterator to loop through lines of a curves.* file each containing
    data fields as defined the in the ecnf-format.txt file, yielding
    pairs (K,ainvs).  Quicker than read_curves as the elliptic curves
    are not constructed.

    If only_one is True, skips curves whose 4th data field is
    *not* 1, hence only yielding one curve per isogeny class.

    """
    count=0
    with open(curvefile) as file:
        for L in file.readlines():
            data = L.split()
            if len(data)!=13:
                print("line {} does not have 13 fields, skipping".format(L))
                continue
            if only_one and data[3]!='1':
                continue
            count +=1
            if ncurves and count>ncurves:
                return
            field_label = data[0]
            K = nf_lookup(field_label)
            conductor_label = data[1]
            iso_label = data[2]
            c_num = data[3]
            conductor_ideal = data[4]
            ainvs = ainvs_from_strings(K, data[6:11])
            yield (field_label,conductor_label,conductor_ideal,iso_label,c_num,ainvs)

def add_trace_hashes(curves_file, isoclass_file, suffix='x', verbose=False):
    r"""One-off function to read a curves file and an isoclass file,
    compute the trace-hashes of one curve in each isogeny class and
    add that to the isoclass file.
    """
    from trace_hash import TraceHash_from_ainvs
    hash_table = {}
    n = 0
    for (field_label,conductor_label,conductor_ideal,iso_label,c_num,ainvs) in read_ai(curves_file, only_one=True):
        label = "-".join([field_label,conductor_label,iso_label])
        hash_table[label+"1"] = ainvs
        n += 1
    print("Finished reading {} curves, now computing hashes".format(n))
    n = 0
    with open(isoclass_file) as isofile, open(isoclass_file+suffix,'w') as new_isofile:
        for L in isofile.readlines():
            label, record = parse_line_label_cols(L)
            if label in hash_table:
                if verbose:
                    print("Processing {}".format(label))
                hash = str(TraceHash_from_ainvs(hash_table[label]))
                new_isofile.write(L[:-1]+" "+hash+"\n")
            else:
                print("label {} in {} has no curve in {}".format(label, isoclass_file, curves_file))
            n+=1
    print("Finished rewriting {} lines to new isoclass file {}".format(n, isoclass_file+"x"))

bad_fields = ['2.2.{}.1'.format(d) for d in [204, 205, 220, 221, 232]]

def fix_models_field(field_type, field_label, verbose=True):
    """Apply fix_model to all curves over one field.  Output new versions
    of curves.<field>, local_data.<field>, mwdata.<field> with ".new"
    suffix.
    """
    from nfscripts import fix_model
    base_dir = os.path.join(ECNF_DIR, field_type)
    data = read_all_field_data(base_dir, field_label, check_cols=False, mwdata_format='new')
    K = nf_lookup(field_label).change_names('w')
    n = 0
    if verbose:
        print("processing curves...")
    for label, record in data.items():
        n += 1
        if verbose and n%1000==0:
            print("{} done...".format(n), end="")
            sys.stdout.flush()
            if n%10000==0:
                print()
        record = fix_model(K, record)
    if(verbose):
        print("...done, writing new files...")

    for ftype in ['curves', 'local_data', 'mwdata']:
        new_file = os.path.join(base_dir, "{}.{}.fixed".format(ftype, field_label))
        with open(new_file, 'w') as outfile:
            n = 0
            for label, record in data.items():
                line = file_line(ftype, record)
                outfile.write(line.rstrip()+"\n")
                n += 1
                if verbose and n%1000==0:
                    print("{} lines output to {}...".format(n, new_file))
    return data

def simplify_ideal_strings_field(field_type, field_label, verbose=True):
    """Convert ideal strings from long form [N,a,alpha] to 1-generator
    (gen) or 2-generatorr (gen1,gen2).  This affects conductor_ideal
    in curves.* and minD, bad_primes in local_data.*.

    At the same time we add the columns 'disc', 'normdisc' to curves.*.

    """
    from nfscripts import simplify_ideal_strings
    base_dir = os.path.join(ECNF_DIR, field_type)
    data = read_all_field_data(base_dir, field_label, check_cols=False, mwdata_format='new')
    K = nf_lookup(field_label)
    if K: # special treatment when 'field_label' is not an actual
          # field label, e.g. for sengun/EC_IQF_egr
        K = K.change_names('w')
    n = 0
    if verbose:
        print("processing ideals...")
    for label, record in data.items():
        n += 1
        if verbose and n%1000==0:
            print("{} done...".format(n), end="")
            sys.stdout.flush()
            if n%10000==0:
                print()
        record = simplify_ideal_strings(K, record)
    if(verbose):
        print("done")
    if verbose:
        print("Writing new files...")

    for ftype in ['curves', 'local_data']: #, 'isoclass', 'mwdata', 'galrep']:
        new_file = os.path.join(base_dir, "{}.{}.new".format(ftype, field_label))
        with open(new_file, 'w') as outfile:
            n = 0
            for label, record in data.items():
                if ftype=='isoclass' and record['number']!=1:
                    continue
                line = file_line(ftype, record)
                outfile.write(line.rstrip()+"\n")
                n += 1
                if verbose and n%1000==0:
                    print("{} lines output to {}...".format(n, new_file))

def convert_curve_file(infilename, outfilename, ncurves=0):
    r"""
    Iterator to loop through lines of a curves.* file each
    containing 13 data fields as defined the in the ecnf-format.txt file,
    writing output file in format needed by Sutherlands galdata program.
    """
    count=0
    pols = {} # keys are fields, values encoded defining polys
    with open(infilename) as file, open(outfilename, mode='w') as outfile:
        for L in file.readlines():
            #stdout.write(L)
            data = L.split()
            if len(data)!=13:
                print("line {} does not have 13 fields, skipping".format(L))
                continue
            count +=1
            if ncurves and count>ncurves:
                outfile.close()
                return

            K = nf_lookup(data[0])
            if not K in pols:
                pols[K] = ",".join([str(c) for c in list(K.defining_polynomial())])
            pol = pols[K]
            label = "-".join(data[:3]) + data[3]
            coeffs = ";".join(data[6:11])
            outfile.write(":".join([label,coeffs,pol])+"\n")

def Q_curve_check(ftypes=all_ftypes, fields=None, certs=False, Detail=1):
    from Qcurves import is_Q_curve
    for ftype in ftypes:
        if Detail:
            print("Checking curves over fields in {}".format(ftype))
        if not fields:
            with open("{}/{}/fields.txt".format(ECNF_DIR,ftype)) as field_file:
                fields = [f[:-1] for f in field_file.readlines()]
        for fname in fields:
            if Detail>1:
                print("Checking curves over field {}".format(fname))
            K = nf_lookup(fname)
            data = {}
            n = 0
            curves_filename = "{}/{}/curves.{}".format(ECNF_DIR,ftype,fname)

            with open(curves_filename) as curves:
                for L in curves.readlines():
                    label, record = parse_curves_line(L)
                    if label:
                        n += 1
                        data[label] = record
            if Detail:
                print("Read {} curves from {}".format(n,curves_filename))
            fcerts = {}
            levels = []
            bads = []
            ngood = 0
            n1 = 0
            for c in data:
                Edata = data[c]
                if Edata['number']!=1:
                    continue
                n1 += 1
                E = curve_from_string(K,Edata['ainvs'])
                lab = Edata['label']
                if Detail>1:
                    print("Running Q-curve test on {}".format(lab))
                res, cert = is_Q_curve(E, certificate=certs, verbose=(Detail>2))
                if res != Edata['q_curve']:
                    print("**********bad result for {}".format(lab))
                    return E, Edata
                    bads.append(lab)
                else:
                    if res:
                        fcerts[lab] = cert
                        if Detail>1:
                            print("yes: certificate = {}".format(cert))
                        if not cert['CM']:
                            N = cert['N']
                            if not N in levels:
                                levels.append(N)
                        if Detail>2:
                            print("{} OK".format(lab))
                    ngood += 1
                if Detail>1 and ngood%100==0:
                    print("{} curves checked OK".format(ngood))
            if bads:
                print("!!!!!!!!! {} discrepancies over {}: {}".format(len(bads), fname, bads))
            if Detail:
                print("Field {}: {} agreements out of {} classes".format(fname, ngood, n1))
                if certs:
                    if fcerts:
                        print("Levels of non-CM Q-curves: {}".format(levels))
                        print("Certificates of Q-curves:")
                        for lab in fcerts:
                            print("{}: {}".format(lab,fcerts[lab]))
                    else:
                        print("No Q-curves")

all_tr_ftypes = ['cubics', 'quartics', 'quintics', 'sextics']

# RQF omitted since we know all curves over those ields are modular

def modularity_check(ftypes=all_tr_ftypes, fields=None, verbose=2):
    for ftype in ftypes:
        if verbose:
            print("Checking curves over fields in {}".format(ftype))
        if fields:
            field_list = fields
        else:
            with open("{}/{}/fields.txt".format(ECNF_DIR,ftype)) as field_file:
                field_list = [f[:-1] for f in field_file.readlines()]

        for fname in field_list:
            if verbose>1:
                print("Checking curves over field {}".format(fname))
            K = nf_lookup(fname)
            n = 0
            curves_filename = "{}/{}/curves.{}".format(ECNF_DIR,ftype,fname)

            for (field_label,conductor_label,conductor_ideal,iso_label,c_num,E) in read_curves(curves_filename):
                n += 1
                lab = "{}-{}-{}{}".format(field_label, conductor_label, iso_label, c_num)
                if verbose>2:
                    print("Checking modularity of {}:".format(lab), end=" ")
                res = isModular(E)
                if res:
                    if verbose>2:
                        print("OK")
                else:
                    print("Modularity check failed for {}".format(lab))
                if n%1000==0 and verbose>1:
                    print("Checked {} curves...".format(n))

def make_mwdata(curves_filename, mwdata_filename, label=None,
                min_cond_norm=None, max_cond_norm=None,
                test_saturation=False, verbose=False):
    r"""Retrieves curves from a curves file
    with conductor norm between given bounds (optional), finds their
    ranks (or bounds) and generators, and outputs an mwdata file.

    If label is given, it should be a short isogeny class label, and
    then only that class will be run.  Otherwise, the minimum and
    maximum conductor may optionally be given.  Otherwise all the
    curves (isogeny classes) in the in put file are processed.

    """
    from mwinfo import get_generators
    from codec import file_line
    with open(mwdata_filename, 'w', 1) as mw_out:
        for cl in read_classes_new(curves_filename):
            short_class_label = "-".join([cl['conductor_label'],cl['iso_label']])
            class_label = "-".join([cl['field_label'],cl['conductor_label'],cl['iso_label']])
            if label:
                if label!=short_class_label:
                    # if verbose:
                    #     print("Skipping {}".format(short_class_label))
                    continue
            NN = cl['conductor_norm']
            if min_cond_norm and NN<min_cond_norm:
                # if verbose:
                #     print("Skipping class as conductor norm < {}".format(min_cond_norm))
                continue
            if max_cond_norm and NN>max_cond_norm:
                if verbose:
                    print("Skipping rest of file as conductor norms >= {} > {}".format(NN,max_cond_norm))
                break

            print("Processing class {}".format(class_label))
            try:
                cl = get_generators(cl, test_saturation=test_saturation, verbose=verbose)
                #mwlines = make_mwdata_lines(cl)
                line = file_line('mwdata', cl)
                if verbose:
                    print(mwlines)
                mw_out.write(mwlines+"\n")
            except RuntimeError as e:
                print("caught RuntimeError: {}".format(e))


def add_analytic_ranks(curves_filename, mwdata_filename, suffix='x', verbose=False):
    r"""Retrieves curves from a curves file and mwdata from the mwdata
     file.  Computes analytic ranks and rewrites the mwdata file
     adding the suffix to its filename.

    This is a one-off since the orginal mwdata file code forgot to
    compute and output analytic ranks.
    """
    from magma import get_magma
    curve_table = {}
    n = 0
    with open(curves_filename) as curves:
        for L in curves.readlines():
            label, record = parse_curves_line(L)
            if label:
                n += 1
                curve_table[label] = record
    print("Read {} curves from {}".format(n,curves_filename))

    ar_table = {}

    def AnalyticRank(Edata):
        nonlocal verbose
        class_label = "-".join([Edata['field_label'],Edata['conductor_label'],Edata['iso_label']])
        if class_label in ar_table:
            return ar_table[class_label]
        if verbose:
            print("Computing analytic rank of {}".format(class_label))
        K = nf_lookup(Edata['field_label'])
        E = curve_from_string(K, Edata['ainvs'])
        magma = get_magma()
        ar = int(magma(E).AnalyticRank())
        ar_table[class_label] = ar
        return ar

    with open(mwdata_filename) as mw_in, open(mwdata_filename+suffix, 'w', 1) as mw_out:
        for L in mw_in.readlines():
            label, record = parse_mwdata_line(L)
            if record['analytic_rank'] is None:
                ar = AnalyticRank(curve_table[label])
                if verbose:
                    print("Updating analytic rank of {} to {}".format(label,ar))
                data = L.split()
                data[6] = str(ar)
                L = " ".join(data)
                if verbose:
                    print("New mwdata line: {}".format(L))
                mw_out.write(L + "\n")
            else:
                mw_out.write(L)

def recompute_real_data(base_dir, field_label, suffix='x', minN=None, maxN=None, one_label=None, max_sat_prime = None, prec=None, Lprec=None, verbose=0):
    r"""
    Reads curves and local data files.
    Computes: analytic rank and L-value using Magma;
     omega (global period);
     heights and regulator;
     analytic Sha (rounded);
    Rewrites mwdata.

    The prec parameter controls the precision to which the heights,
    regulator and global period is computed.  It is bit precision.
    The Lprec parameter (also bit precision) controls the precision
    used for the L-value in Magma: note that the running time increases
    rapidly!

        # prec (bits) magma_prec (decimal)
        15-18 5
        19-21 6
        22-24 7
        25-28 8
        29-31 9
        32-34 10
        35-38 11
        39-41 12
        42-44 13
        45-48 14
        49-51 15
        52-54 16
        55-58 17
        59-61 18
        62-64 19
        65-68 20
       128-131 39

    The real work is done in the function extend_mwdata_one() from
    nfscripts.py.  This function is similar to extend_mwdata().

    """
    from magma import get_magma
    from nfscripts import extend_mwdata_one

    data = read_all_field_data(base_dir, field_label, check_cols=False)
    classdata = {} # will hold isogeny-invariant values keyed by class label
    Kfactors = {} # BSD factor depending only on the field K

    if one_label:
        mwoutfile = base_dir+'/mwdata.'+suffix+"."+one_label
    else:
        mwoutfile = base_dir+'/mwdata.'+field_label+suffix

    with open(mwoutfile, 'w', 1) as mwdata:
        for label, Edata in data.items():
            class_label = Edata['class_label']
            if one_label and one_label!=class_label:
                continue
            N = Edata['conductor_norm']
            if (minN and N<minN) or (maxN and N>maxN):
                continue
            print("Processing {}".format(label))
            if not class_label in classdata:
                # then this is a new isogeny class, so we'll use magma in extend_mwdata_one
                magma = get_magma()

            Edata = extend_mwdata_one(Edata, classdata, Kfactors, magma,
                                      max_sat_prime = max_sat_prime,
                                      prec=prec, Lprec=Lprec, verbose=verbose)
            line = file_line('mwdata', Edata)
            if verbose>1:
                print("New mwdata line: {}".format(line))
            mwdata.write(line + "\n")

def write_data_files(data, file_types=all_file_types, field_type=None, field_label='test', base_dir=ECNF_DIR, append=False, suffix='part'):
    """
    data is a dict whose values are curve records.

    Outputs to files basedir/field_type/<ft>.field_name, either appending or overwriting.
    """
    from schemas import column_names
    mode = 'a' if append else 'w'
    for ft in file_types:
        new_file = os.path.join(base_dir, field_type, f"{ft}.{field_label}.{suffix}")
        cols = column_names[ft]
        with open(new_file, mode) as outfile:
            n = 0
            for label, record in data.items():
                if ft=='isoclass' and record['number']!=1:
                    continue
                if not all(col in record for col in cols):
                    print("Incomplete record for {}".format(label))
                    print("Expected keys: {}".format(cols))
                    print("Missing keys: {}".format([col for col in cols if col not in record]))
                line = file_line(ft, record)
                outfile.write(line.rstrip()+"\n")
                n += 1
                if n%1000==0:
                    print("{} lines output to {} so far...".format(n, new_file))
        print("{} lines output to {} so far...".format(n, new_file))

def make_all_data_files(raw_curves, file_types=all_file_types,
                        field_type=None, field_label='test', base_dir=ECNF_DIR, verbose=0,
                        prec=None):
    """raw_curves is a generator yielding short 'raw' curve dicts with
    fields (as in read_curves_magma()):

    field_label, conductor_norm, conductor_label, conductor_ideal, iso_label, ainvs

    This computes isogeny classes and all curve data for all curves in
    each class, writing the result to files
    <base_dir>/<field_type>/<ft>.<field_label> for each file typ ft
    (default: all, i.e. curves, isoclass, local_data, mwdata, galrep)

    The full data is also returned.

    prec controls the bit precision of heights, special L-value, etc.
    If None (the default) is uses standard 53-bit precision.

    """
    from nfscripts import make_isogeny_classes
    data = make_isogeny_classes(raw_curves, verbose=verbose, prec=prec)
    write_data_files(data, file_types, field_type, field_label, base_dir)
    return data

def make_all_data_files1(raw_curves, file_types=all_file_types,
                        field_type=None, field_label='test', base_dir=ECNF_DIR, verbose=0,
                        prec=None):
    """raw_curves is a generator yielding short 'raw' curve dicts with
    fields (as in read_curves_magma()):

    field_label, conductor_norm, conductor_label, conductor_ideal, iso_label, ainvs

    This computes isogeny classes and all curve data for all curves in
    each class, writing the result to files
    <base_dir>/<field_type>/<ft>.<field_label> for each file typ ft
    (default: all, i.e. curves, isoclass, local_data, mwdata, galrep)

    prec controls the bit precision of heights, special L-value, etc.
    If None (the default) is uses standard 53-bit precision.

    Unlike make_all_data_files() which calls make_isogeny_classes just
    once and only outputs right at the end, this one callas
    make_isogeny_class() for each raw curve and outputs as it goes
    along.  This guards against sage/magma/pari crashing in the middle
    of a long run.

    """
    from nfscripts import make_isogeny_class
    for curve in raw_curves:
        label = "{}-{}-{}".format(curve['field_label'],curve['conductor_label'],curve['iso_label'])
        print("working on class {}".format(label))
        data = make_isogeny_class(curve, verbose=verbose, prec=prec)
        write_data_files(data, file_types, field_type, field_label, base_dir, append=True)
        print("output for class {} complete".format(label))
        print("====================================")
