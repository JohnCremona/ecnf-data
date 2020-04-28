# Functions for reading and writing data files

import re
from sage.all import ZZ, QQ, latex, Set
from sage.databases.cremona import cremona_to_lmfdb
from codec import convert_conductor_label, curve_from_string, curve_from_strings, ainvs_from_strings, convert_ideal_label, local_data_to_string, ideal_to_string, NFelt, curves_data_to_string, curve_from_data, local_data_from_string
from fields import nf_lookup
from os import getenv
#from sys import stdout

# Functions to parse a single line from one of the files curves.*, isoclass.*, mwdata.*, local_data.*, galdata.*
#
# In each case the function returns a full label and a dict whose kets are exactly the relevant table columns
#

# The first 4 columns in curves.*, isoclass.*, mwdata.*, local_data.*
# are the same and define the full label, so we factor out this part.
#

def split_field_label(field_label):
    r"""Return (degree, signature, abs(disc), index) from an LMFDB field label
    """
    d, r, a, i = field_label.split(".")
    d = int(d)
    r = int(r)
    a = ZZ(a)
    i = int(i)
    s = [r,(d-r)//2]
    return (d,s,a,i)

def numerify_iso_label(lab):
    r"""Return the numerical equivalent of an isogeny class letter-code.

    Normally this is the integer>=0 which lab represents in base 26,
    with digits a=0,...,z=25, but if lab starts with 'CM' and the rest
    represents the integer n as above, the numerical version is
    -n-1. (Not -n as then 0 would be ambiguous).  This variant is only
    relevant (currently) over imaginary quadratic fields.  The label
    may be in upper case (currently only over 3.1.23.1), but that
    should not be used over any field for whcih a CM label is
    possible.
    """
    from sage.databases.cremona import class_to_int
    if 'CM' in lab:
        return -1 - class_to_int(lab[2:])
    else:
        return class_to_int(lab.lower())

def split_galois_image_code(s):
    """Each code starts with a prime (1-3 digits but we allow for more)
    followed by an image code for that prime.  This function returns
    two substrings, the prefix number and the rest.
    """
    p = re.findall(r'\d+', s)[0]
    return p, s[len(p):]


def parse_line_label_cols(L):
    r"""
    Parse the first 4 columns of one line from a curves/isoclass/mwdata/local_data file
    """
    data = L.split()
    record = {}
    record['field_label'] = field_label = data[0]
    degree, signature, abs_disc, _ = split_field_label(field_label)
    record['degree'] = degree
    record['signature'] = signature
    record['abs_disc'] = abs_disc

    conductor_label = data[1]
    iso_label = data[2]
    number = int(data[3])
    short_label = "{}-{}{}".format(conductor_label, iso_label, number)
    label = "{}-{}".format(field_label, short_label)
    short_class_label = "{}-{}".format(conductor_label, iso_label)

    record['conductor_label'] = conductor_label
    record['iso_label'] = iso_label
    record['iso_nlabel'] = numerify_iso_label(iso_label)
    record['number'] = number
    record['short_label'] = short_label
    record['label'] = label
    record['short_class_label'] = short_class_label
    record['class_label'] = "{}-{}".format(field_label, short_class_label)
    return label, record

def parse_curves_line(L):
    r"""
    Parse one line from a curves file
    """
    data = L.split()
    if len(data)!=12:
        print("curves line {} does not have 12 fields, skipping".format(L))
        return
    label, record = parse_line_label_cols(L)

    record['conductor_ideal'] = data[4]
    record['conductor_norm'] = ZZ(data[5])

    record['ainvs'] = data[6]
    record['jinv'] = data[7]
    record['equation'] = data[8]
    record['cm'] = ZZ(data[9]) if data[9]!='?' else '?'
    bc = data[10][1:-1]
    record['base_change'] = [str(lab) for lab in bc.split(",")] if bc else []
    record['q_curve'] = (data[1]==1)

    return label, record

def parse_isoclass_line(L):
    r"""
    Parse one line from an isoclass file
    """
    data = L.split()
    if len(data)!=6:
        print("isoclass line {} does not have 6 fields, skipping".format(L))
        return
    label, record = parse_line_label_cols(L)

    record['isogeny_matrix'] = mat = [[int(a) for a in r.split(",")]
                                      for r in data[4][2:-2].split("],[")]
    record['class_size'] = len(mat)
    record['class_deg'] = max(max(r) for r in mat)
    record['all_iso_degs'] = dict([[n+1,sorted(list(set(row)))] for n,row in enumerate(mat)])
    record['trace_hash'] = ZZ(data[5])

    # NB Every curve in the class has the same 'isogeny_matrix',
    # 'class_size', 'class_deg', and the for the i'th curve in the
    # class (for i=1,2,3,...) its 'isogeny_degrees' column is
    # all_iso_degs[i].

    return label, record

def parse_local_data_line(L):
    r"""
    Parse one line from a local_data file
    """
    data = L.split()
    ncols = len(data)
    if not ncols in [6,7]:
        print("local_data line {} does not have 6 or 7 fields, skipping".format(L))
        return
    label, record = parse_line_label_cols(L)

    record['local_data'] = [] if (ncols==6) else local_data_from_string(data[4])
    # The non_min_p column is a list of strings
    # e.g. ['[N1,a1,alpha1]', '[N2,a2,alpha2]'] while the string in
    # the file will contain [[N1,a1,alpha1],[N2,a2,alpha2]].
    # Currently the list has 0 or 1 entries but we do not want to rely
    # on this.
    nmp = data[-2]
    record['non_min_p'] = [] if nmp == '[]' else ['['+id+']' for id in nmp[2:-2].split("],[")]
    record['minD'] = data[-1]

    return label, record

def parse_mwdata_line(L):
    r"""
    Parse one line from an mwdata file
    """
    data = L.split()
    if len(data)!=14:
        print("mwdata line {} does not have 14 fields, skipping".format(L))
        return
    label, record = parse_line_label_cols(L)

    r = data[4]
    record['rank'] = None if r=='?' else int(r)
    r = data[5]
    record['rank_bounds'] = '?' if r=='?' else [int(rb) for rb in r[1:-1].split(",")]
    r = data[6]
    record['analytic_rank'] = None if r=='?' else int(r)
    record['ngens'] = int(data[7])
    gens = data[8]
    record['gens'] = [] if gens == '[]' else gens.replace("[[[","[[").replace("]]]","]]").replace("]],[[","]];[[").split(";")
    record['heights'] = data[9]
    record['reg'] = data[10]
    record['torsion_order'] = int(data[11])
    ts = data[12]
    record['torsion_structure'] = [] if ts=='[]' else [int(t) for t in ts[1:-1].split(",")]
    tgens = data[13]
    record['torsion_gens'] = [] if tgens == '[]' else tgens.replace("[[[","[[").replace("]]]","]]").replace("]],[[","]];[[").split(";")

    return label, record

def parse_galrep_line(L):
    r"""
    Parse one line from a galrep file
    """
    data = L.split()
    label = data[0]
    galois_images = data[1:]
    pr = [ int(split_galois_image_code(s)[0]) for s in galois_images]
    record = {'label': label,
              'non-surjective_primes': pr,
              'galois_images': galois_images,
    }

    return label, record

# Python types of the columns in ec_nfcurves:

from six import text_type
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
                  'ngens': int_type,
                  'rank': int_type,
                  'rank_bounds': list_type, # 2 ints
                  'analytic_rank': int_type,
                  'torsion_order': int_type,
                  'torsion_structure': list_type, # 0,1,2 ints
                  'gens': list_type, # of strings
                  'torsion_gens': list_type, # of strings
                  'isogeny_matrix': list_type, # of lists of ints
                  'isogeny_degrees': list_type, # of ints
                  'class_deg': int_type,
                  'non-surjective_primes': list_type, # of ints
                  'galois_images': list_type, # of strings
                  'equation': str_type,
                  'local_data': list_type, # of dicts
                  'non_min_p': list_type, # of strings
                  'minD': str_type,
                  'heights': list_type, # of floats
                  'reg': float_type, # or int(1)
                  'q_curve': bool_type,
                  'base_change': list_type, # of strings
                  'trace_hash': hash_type
}

ec_nfcurves_columns = Set(['field_label', 'degree', 'signature',
                           'abs_disc', 'label', 'short_label',
                           'class_label', 'short_class_label',
                           'conductor_label', 'iso_label',
                           'iso_nlabel', 'conductor_ideal',
                           'conductor_norm', 'number', 'ainvs',
                           'jinv', 'equation', 'cm', 'base_change',
                           'q_curve', 'class_size', 'isogeny_matrix',
                           'class_deg', 'isogeny_degrees',
                           'trace_hash', 'local_data', 'non_min_p',
                           'minD', 'rank', 'rank_bounds',
                           'analytic_rank', 'ngens', 'gens',
                           'heights', 'reg', 'torsion_order',
                           'torsion_structure', 'torsion_gens',
                           'galois_images', 'non-surjective_primes'])

assert ec_nfcurves_columns==Set(keys_and_types.keys())

def read_all_field_data(base_dir, field_label, check_cols=True):
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
    curves_filename = '{}/curves.{}'.format(base_dir,field_label)
    isoclass_filename = '{}/isoclass.{}'.format(base_dir,field_label)
    local_data_filename = '{}/local_data.{}'.format(base_dir,field_label)
    mwdata_filename = '{}/mwdata.{}'.format(base_dir,field_label)
    galrep_filename = '{}/galrep.{}'.format(base_dir,field_label)

    all_data = {}
    n = 0

    with open(curves_filename) as curves:
        for L in curves.readlines():
            label, record = parse_curves_line(L)
            if label:
                n += 1
                all_data[label] = record
    print("Read {} curves from {}".format(n,curves_filename))
    n = 0

    with open(isoclass_filename) as isoclass:
        for L in isoclass.readlines():
            label, record = parse_isoclass_line(L)
            if label:
                n += 1
                nc = record['class_size']
                all_iso_degs = record.pop('all_iso_degs')
                label = label[:-1] # delete final '1'
                for ic in range(1,nc+1):
                    record['isogeny_degrees'] = all_iso_degs[ic]
                    all_data[label+str(ic)].update(record)
    print("Read {} classes from {}".format(n,isoclass_filename))
    n = 0

    with open(local_data_filename) as local_data:
        for L in local_data.readlines():
            label, record = parse_local_data_line(L)
            if label:
                n += 1
                all_data[label].update(record)
    print("Read {} local_data records from {}".format(n,local_data_filename))
    n = 0

    with open(mwdata_filename) as mwdata:
        for L in mwdata.readlines():
            label, record = parse_mwdata_line(L)
            if label:
                n += 1
                all_data[label].update(record)
    print("Read {} mwdata records from {}".format(n,mwdata_filename))
    n = 0

    with open(galrep_filename) as galrep:
        for L in galrep.readlines():
            label, record = parse_galrep_line(L)
            if label:
                n += 1
                all_data[label].update(record)
    print("Read {} galrep records from {}".format(n,galrep_filename))

    if check_cols:
        for label in all_data:
            cols = Set(all_data[label])
            if cols != ec_nfcurves_columns:
                print("Wrong key set for {}".format(label))
                diff = cols - ec_nfcurves_columns
                if diff:
                    print("data has extra keys {}".format(diff))
                diff = ec_nfcurves_columns - cols
                if diff:
                    print("data is missing keys {}".format(diff))
    return all_data

def read_curve_file(infile):
    r"""
    Read a curves file, each line containing 13 data fields as defined
    the in the ecnf-format.txt file. Output is a list of dicts, one
    per curve, in the same order as input.
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
                     'N_label': data[1],
                     'iso_label': data[2],
                     'c_num': data[3],
                     'N_def': data[4],
                     'N_norm': data[5],
                     'ainvs': data[6:11],
                     'cm_flag': data[11],
                     'q_curve_flag': data[12]
            }
            yield curve
    return

def read_curvedata_file(infile):
    r"""
    Read a curvedata file, each line containing 9+ data fields as defined
    the in the ecnf-format.txt file. Output is a list of dicts, one
    per curve, in the same order as input.
    """
    index = 0
    with open(infile) as file:
        for L in file.readlines():
            if L[0]=='#': # allow for comment lines
                continue
            data = L.split()
            ngens = int(data[7])
            if len(data) != 9+ngens:
                print("line {} does not have 9+{} fields, skipping".format(L, ngens))
                continue
            index += 1
            curve = {'index': index,
                     'field_label': data[0],
                     'N_label': data[1],
                     'iso_label': data[2],
                     'c_num': data[3],
                     'rank': data[4],
                     'rank_bounds': data[5],
                     'an_rank': data[6],
                     'ngens': data[7],
                     'gens': data[8:8+ngens],
                     'sha_an': data[8+ngens]
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
                     'N_label': data[1],
                     'iso_label': data[2],
                     'c_num': data[3],
                     'isomat': data[4]
            }
            yield curve
    return

def write_curve_file(curves, outfile):
    r"""
    Write a curves file, each line containing 13 data fields as defined
    the in the ecnf-format.txt file. Input is a list of dicts.
    """
    with open(outfile, 'w') as out:
        for c in curves:
            line = " ".join([c['field_label'],
                             c['N_label'],
                             c['iso_label'],
                             c['c_num'],
                             c['N_def'],
                             c['N_norm']] + c['ainvs'] + [c['cm_flag'],
                                                          c['q_curve_flag']])
            out.write(line+"\n")

def write_curvedata_file(curves, outfile):
    r"""
    Write a curvedata file, each line containing 9+ngens data fields as defined
    the in the ecnf-format.txt file. Input is a list of dicts.
    """
    with open(outfile, 'w') as out:
        for c in curves:
            line = " ".join([c['field_label'],
                             c['N_label'],
                             c['iso_label'],
                             c['c_num'],
                             c['rank'],
                             c['rank_bounds'],
                             c['an_rank'],
                             c['ngens']] + c['gens'] + [c['sha_an']])
            out.write(line+"\n")

def write_isoclass_file(curves, outfile):
    r"""
    Write an isoclass file, each line containing 5 data fields as defined
    the in the ecnf-format.txt file. Input is a list of dicts.
    """
    with open(outfile, 'w') as out:
        for c in curves:
            line = " ".join([c['field_label'],
                             c['N_label'],
                             c['iso_label'],
                             c['c_num'],
            c['isomat']])
            out.write(line+"\n")

def rewrite_curve_file(infile, outfile, verbose=True):
    """
    Convert ideal labels for IQFs fro old-atyle N.c.d to new N.i LMFDB
    standard.  Could also be used over other fields to standardise
    labels into the canonical LMFDB ordering of ideals of the same
    norm.

    Can be used for curves* files,  curve_data* and isoclass* files.
    """
    if 'curves' in infile:
        read_file = read_curve_file
        write_file = write_curve_file
    elif 'curve_data' in infile:
        read_file = read_curvedata_file
        write_file = write_curvedata_file
    elif 'isoclass' in infile:
        read_file = read_isoclass_file
        write_file = write_isoclass_file
    else:
        print("Invalid filename {}: should be a curves or curve_data or isoclass file")
        return

    def convert(c):
        field_label = c['field_label']
        cond_label =  c['N_label']
        new_cond_label = convert_conductor_label(field_label, cond_label)
        if verbose:
            print("Conductor label {} converted to {}".format(cond_label, new_cond_label))
        c['N_label'] = new_cond_label
        return c

    def converter():
        for c in read_file(infile):
            yield convert(c)
        return

    write_file(converter(), outfile)

def read_curves(infile, only_one=False, ncurves=0):
    r"""
    Iterator to loop through lines of a curves.* file each
    containing 13 data fields as defined the in the ecnf-format.txt file,
    yielding its curves as EllipticCurve objects.

    If only_one is True, skips curves whose 4th data field is
    *not* 1, hence only yielding one curve per isogeny class.
    """
    count=0
    with open(infile) as file:
        for L in file.readlines():
            #stdout.write(L)
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
            N_label = data[1]
            iso_label = data[2]
            c_num = data[3]
            N_def = data[4]
            E = curve_from_strings(K, data[6:11])
            yield (field_label,N_label,N_def,iso_label,c_num,E)

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

def read_classes(infile):
    r"""
    Iterator to loop through lines of a curves.* file each
    containing 13 data fields as defined the in the ecnf-format.txt file,
    yielding complete isogeny classes.  These are dicts with keys

    field_label, N_label, N_def, N_norm, iso_label, curves

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
            N_label = data[1]
            iso_label = data[2]
            #c_num = data[3]
            N_def = data[4]
            N_norm = int(data[5])
            E = curve_from_strings(K, data[6:11])
            this_class_id = "-".join(data[:3])
            if this_class_id == prev_class_id:
                this_class['curves'].append(E)
            else:
                if this_class: # else we have just started
                    yield this_class
                this_class = {}
                this_class['field_label'] = field_label
                this_class['N_label'] = N_label
                this_class['N_def'] = N_def
                this_class['N_norm'] = N_norm
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

    field_label, N_label, N_def, N_norm, iso_label, curves

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
            N_label = record['conductor_label']
            iso_label = record['iso_label']
            N_def = record['conductor_ideal']
            N_norm = record['conductor_norm']
            E = curve_from_string(K, record['ainvs'])
            this_class_id = record['class_label']
            if this_class_id == prev_class_id:
                this_class['curves'].append(E)
            else:
                if this_class: # else we have just started
                    yield this_class
                this_class = {}
                this_class['field_label'] = field_label
                this_class['N_label'] = N_label
                this_class['N_def'] = N_def
                this_class['N_norm'] = N_norm
                this_class['iso_label'] = iso_label
                this_class['curves'] = [E]
                prev_class_id = this_class_id
    # final class:
    yield this_class

HOME = getenv("HOME")
bianchi_data_dir = HOME + "/bianchi-data"

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
        #if verbose:
        #    print("raw input: %s" % L)
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
        nf['aq'] = [int(e) for e in ALs[1:-1].split(",")]
        nf['ap'] = [int(e) for e in aplist[1:-1].split(",")]
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
    out = open(outfile, mode='w')
    for L in open(bianchi_data_dir + "/ideals/" + infile).readlines():
        field, ideal = L.split()
        label = convert_ideal_label(nf_lookup(field),ideal)
        out.write(' '.join([field, ideal, label])+'\n')
    out.close()

def make_local_data_file(curves_filename, ld_filename, verbose=False):
    r"""Create a local_data file from a curves file.  This will not be
    needed once we create the local_data file at the same time as the
    curves files.
    """
    from nfscripts import local_data
    with open(ld_filename, 'w', 1) as ldfile:
        for  (field_label,N_label,N_def,iso_label,c_num,E) in read_curves(curves_filename):
            if verbose:
                print("Processing {}".format("-".join([field_label,N_label,iso_label])+c_num))
                #print("E = {}".format(E.ainvs()))
            Eld, nonminP, minD = local_data(E)
            #print("local data: {}".format(Eld))
            Eld = local_data_to_string(Eld)
            #print("local data encoded: {}".format(Eld))
            nonminP = str([ideal_to_string(P) for P in nonminP]).replace(" ","")
            minD = ideal_to_string(minD)
            line = " ".join([field_label,N_label,iso_label,c_num,Eld, nonminP, minD])
            ldfile.write(line + "\n")

def extend_curves_file(infilename, outfilename, verbose=False):
    r"""One-off function to extend a curves file from the old format (13
    columns, with the ai taking 5 columns) to the new (12 columns,
    just one for the ai and three extras), the extra ones being
    'jinv', 'equation', 'base_change'

    In addition to adding the three columns we ensure that the
    is_Q_curve column is not '?'.
    """
    nmax = -1 # make this positive for a test run which stops after this many lines
    n = 0
    with open(outfilename, 'w', 1) as out:
        for curve in read_curve_file(infilename):
            n += 1
            label = "-".join([curve['field_label'], curve['N_label'], curve['iso_label']]) + curve['c_num']
            if verbose:
                print("Processing {}".format(label))
                line = curves_data_to_string(curve, old_style=True)
                print("Old line:\n{}".format(line))
            E = curve_from_data(curve)
            curve['ainvs'] = ";".join(curve['ainvs'])
            curve['jinv'] = NFelt(E.j_invariant())
            curve['equation'] = str(latex(E)).replace(" ","") # no "\(", "\)"
            if curve['q_curve_flag'] == '?':
                from nfscripts import is_Q_curve
                qc = is_Q_curve(E)
                print("Filling in Q-curve flag for curve {} to {}".format(label, qc))
                curve['q_curve_flag'] = int(qc)
            EQlist = E.descend_to(QQ)
            if EQlist:
                bc = [cremona_to_lmfdb(EQ.label()) for EQ in EQlist]
                if verbose:
                    print("{} is base change of {}".format(label, bc))
            else:
                bc = []
            curve['base_change'] = "[" + ",".join(bc) + "]"
            if False:
                print("New curve struct: {}".format(curve))
            line = curves_data_to_string(curve)
            if verbose:
                print("New line:\n{}".format(line))
            out.write(line+"\n")
            if n==nmax:
                break
            if n%1000==0:
                print("{} curves processed".format(n))

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
            #stdout.write(L)
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
            N_label = data[1]
            iso_label = data[2]
            c_num = data[3]
            N_def = data[4]
            ainvs = ainvs_from_strings(K, data[6:11])
            yield (field_label,N_label,N_def,iso_label,c_num,ainvs)

def add_trace_hashes(curves_file, isoclass_file, suffix='x', verbose=False):
    r"""One-off function to read a curves file and an isoclass file,
    compute the trace-hashes of one curve in each isogeny class and
    add that to the isoclass file.
    """
    from trace_hash import TraceHash_from_ainvs
    hash_table = {}
    n = 0
    for (field_label,N_label,N_def,iso_label,c_num,ainvs) in read_ai(curves_file, only_one=True):
        label = "-".join([field_label,N_label,iso_label])
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

