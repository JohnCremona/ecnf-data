# Functions for reading and writing data files

import re
from sage.all import ZZ, QQ, latex, Set, Magma, RealField, RR, Infinity
from sage.databases.cremona import cremona_to_lmfdb
from codec import convert_conductor_label, curve_from_string, curve_from_strings, ainvs_from_strings, convert_ideal_label, local_data_to_string, ideal_to_string, NFelt, curves_data_to_string, curve_from_data, local_data_from_string, encode_int_list, decode_int_list, decode_points_one2many, encode_points_many2one, parse_point, encode_points
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

    ldstring = "" if (ncols==6) else data[4]
    ld, ldx = local_data_from_string(ldstring)
    record['local_data'] = ld
    record.update(ldx) # fields 'badp', 'nbadp', 'ss', 'tamprod'

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
    record['torsion_order'] = nt = int(data[11])
    ts = data[12]
    record['torsion_structure'] = [] if ts=='[]' else [int(t) for t in ts[1:-1].split(",")]
    record['torsion_primes'] = ZZ(nt).prime_divisors()
    tgens = data[13]
    record['torsion_gens'] = [] if tgens == '[]' else tgens.replace("[[[","[[").replace("]]]","]]").replace("]],[[","]];[[").split(";")

    return label, record

def parse_new_mwdata_line(L):
    r"""
    Parse one line from an mwdata file (with extra columns omega, lvalue, sha)
    """
    data = L.split()
    if len(data)!=17:
        print("mwdata line {} does not have 14 fields, skipping".format(L))
        return
    label, record = parse_line_label_cols(L)

    def decode_col(col, decoder): # use for columns which may have '?'
        return None if col=='?' else decoder(col)

    record['rank']              = decode_col(data[4], int)
    record['rank_bounds']       = decode_col(data[5], decode_int_list)
    record['analytic_rank']     = decode_col(data[6], int)
    record['ngens']             = int(data[7])
    record['gens']              = decode_points_one2many(data[8])
    record['heights']           = data[9]
    record['reg']               = decode_col(data[10], RR)
    record['torsion_order']     = nt = int(data[11])
    record['torsion_primes']    = ZZ(nt).prime_divisors()
    record['torsion_structure'] = decode_int_list(data[12])
    record['torsion_gens']      = decode_points_one2many(data[13])
    record['omega']             = RR(data[14])
    record['Lvalue']            = RR(data[15])
    record['sha']               = decode_col(data[16], int)

    return label, record

def make_line_label_cols(Edata):
    r"""
    Form string containing the 4 label columns from a curve record
    """
    return " ".join([Edata['field_label'], Edata['conductor_label'], Edata['iso_label'], str(Edata['number'])])

def make_mwdata_line(Edata):
    r"""
    Form one line of an mwdata file from a curve record
    """
    label_cols = make_line_label_cols(Edata)

    rank = str(Edata['rank']) if Edata['rank']!=None else '?'
    reg  = str(Edata['reg'])  if Edata['reg']  else '?'
    sha  = str(Edata['sha'])  if Edata['sha']  else '?'
    rbds = str(Edata['rank_bounds']).replace(" ","")
    gens = encode_points_many2one(Edata['gens'])
    tors = encode_int_list(Edata['torsion_structure'])
    torgens = encode_points_many2one(Edata['torsion_gens'])

    fields = [label_cols,
              rank,
              rbds,
              str(Edata['analytic_rank']),
              str(Edata['ngens']),
              gens,
              Edata['heights'],
              reg,
              str(Edata['torsion_order']),
              tors,
              torgens,
              str(Edata['omega']),
              str(Edata['Lvalue']),
              sha]

    return " ".join(fields)

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

def read_curves_new(infile, only_one=False, ncurves=0):
    r"""
    Iterator to loop through lines of a curves.* file each
    containing 12 data fields as defined the in the ecnf-format.txt file,
    yielding its curves as EllipticCurve objects.

    If only_one is True, skips curves whose 4th data field is
    *not* 1, hence only yielding one curve per isogeny class.
    """
    count=0
    with open(infile) as file:
        for L in file.readlines():
            #stdout.write(L)
            data = L.split()
            if len(data)!=12:
                print("line {} does not have 12 fields, skipping".format(L))
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
            E = curve_from_string(K, data[6])
            yield (field_label,N_label,N_def,iso_label,c_num,E)

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

    The prec parameter only affects the precision to which the global
    period is computed.
    """
    from nfscripts import local_data
    with open(ld_filename, 'w', 1) as ldfile:
        for  (field_label,N_label,N_def,iso_label,c_num,E) in read_curves(curves_filename):
            if verbose:
                print("Processing {}".format("-".join([field_label,N_label,iso_label])+c_num))
            Eld, nonminP, minD = local_data(E)
            Eld = local_data_to_string(Eld)
            nonminP = str([ideal_to_string(P) for P in nonminP]).replace(" ","")
            minD = ideal_to_string(minD)
            line = " ".join([field_label,N_label,iso_label,c_num,Eld, nonminP, minD])
            ldfile.write(line + "\n")

def extend_mwdata(base_dir, field_label, suffix='x', minN=None, maxN=None, one_label=None, max_sat_prime = Infinity, prec=None, verbose=False):
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
    from nfscripts import global_period
    data = read_all_field_data(base_dir, field_label, check_cols=False)
    classdata = {} # will hold isogeny-invariant values keyed by class label
    sat_needed = {} # will hold True/False keyed by class label if curve #1's gens were not saturated

    if prec is None:  # Magma's precision variable is decimal, 53 bits is 16 digits
        RR = RealField()
        prec = RR.precision()
        magma_prec = 16
    else:
        RR  = RealField(prec)
        # log(2)/log(10) =  0.301029995663981
        magma_prec = (prec*0.301029995663981).round()

    Kfactors = {} # BSD factor depending only on the field K
    nmag = 0 # count the number of times we use a Magma instance, and restart every 100
    magma = Magma()
    if one_label:
        mwoutfile = base_dir+'/mwdata.'+field_label+suffix+"."+one_label
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
                #print(Edata)
            K = nf_lookup(Edata['field_label'])
            # We need to construct every E as a Sage EllipticCurve in
            # order to compute omega, but we only need construct it as
            # a Magma curve once per isogeny class.
            E = curve_from_string(K,Edata['ainvs'])

            # find analytic rank and L-value:

            class_label = Edata['class_label']
            if not class_label in classdata:
                nmag += 1
                if nmag%100==0:
                    magma.quit(verbose=verbose)
                    magma = Magma()
                    nmag = 0
                mE = magma(E)
                if verbose:
                    print("Calling Magma's AnalyticRank()")
                ar, lval = mE.AnalyticRank(Precision=magma_prec, nvals=2)
                lval = RR(lval)
                ar = int(ar)
                classdata[class_label] = (ar,lval)
            else:
                ar, lval = classdata[class_label]
            Edata['analytic_rank'] = ar
            Edata['Lvalue'] = lval
            if verbose:
                print("analytic rank = {}\nL-value = {}".format(ar,lval))

            # recompute regulator.  Original heights were computed
            # before fixing Sage's height function precision issues
            # properly.

            gens = [E(parse_point(K,P)) for P in Edata['gens']]
            if verbose:
                print("gens = {}".format(gens))

            if max_sat_prime and ngens and (Edata['number']==1 or sat_needed['class_label']):
                if max_sat_prime==Infinity:
                    new_gens, index, new_reg = E.saturation(gens, verbose=verbose)
                else:
                    new_gens, index, new_reg = E.saturation(gens, max_prime=max_sat_prime)
                if index>1:
                    print("Original gens were not saturated, index = {} (using max_prime {})".format(index,max_sat_prime))
                    gens = new_gens
                    if Edata['number']==1:
                        sat_needed['class_label'] = True
                else:
                    if verbose:
                        print("gens are saturated at primes up to {}".format(max_sat_prime))
                    if Edata['number']==1:
                        sat_needed['class_label'] = False

            heights = [P.height(precision=prec) for P in gens]
            Edata['heights'] = str(heights).replace(" ","")
            if verbose:
                print("heights = {}".format(heights))

            reg = E.regulator_of_points(gens, precision=prec)
            Edata['reg'] = str(reg) if ar else '1'
            if verbose:
                print("regulator (of known points) = {}".format(reg))

            # allow for the scaling in the Neron-Tate height
            # pairing: for BSD we need non-normalised heights and
            # normalization divides every height by K.degree(), so
            # the regulator we need has to be multiplied by
            # K.degree()**rank.
            if len(gens) == ar:
                NTreg = reg * K.absolute_degree()**ar
                if verbose:
                    print("Neron-Tate regulator = {}".format(NTreg))
            else:
                NTreg = None

            # compute omega

            # find scaling factor in case we don't have a global minimal model
            minDnorm = ZZ(Edata['minD'][1:].split(",")[0]).abs()
            modelDnorm = E.discriminant().norm().abs()
            fac = (modelDnorm/minDnorm).nth_root(12) # will be exact
            if fac!=1 and verbose:
                print("Not a global minimal model")
                print("Scaling factor = {}".format(fac))
            Edata['omega'] = omega = global_period(E, fac, prec=prec)
            if verbose:
                print("omega = {}".format(omega))

            T = E.torsion_subgroup()
            nt = T.order()
            if nt != Edata['torsion_order']:
                print("{}: torsion order is {}, not {} as on file; updating data".format(label,nt,Edata['torsion_order']))
                Edata['torsion_order'] = nt
                Edata['torsion_structure'] = list(T.invariants())
                Edata['torsion_gens'] = encode_points([P.element() for P in T.gens()])
            if verbose:
                print("Torsion order = {} (checked)".format(nt))

            tamagawa_product = Edata['tamprod']
            if verbose:
                print("Tamagawa product = {}".format(tamagawa_product))

            if NTreg:
                Rsha = lval * nt**2  / (NTreg * tamagawa_product * omega)

                if not K in Kfactors:
                    Kfactors[K] = RR(K.discriminant().abs()).sqrt() / 2**(K.signature()[1])
                if verbose:
                    print("Field factor = {}".format(Kfactors[K]))

                Rsha *= Kfactors[K]
                Edata['sha'] = sha = Rsha.round()
                if verbose:
                    print("Approximate analytic Sha = {}, rounds to {}".format(Rsha, sha))
                if sha==0 or (sha-Rsha).abs()>0.0001 or not ZZ(sha).is_square():
                    if not verbose:
                        print("Approximate analytic Sha = {}, rounds to {}".format(Rsha, sha))
                    print("****************************Not good! 0 or non-square or not close to a positive integer!")

                # gens = [E(parse_point(K,P)) for P in Edata['gens']]
                # magma = Magma()
                # Rsha_magma = RR(magma(E).ConjecturalSha(magma(gens)))
                # magma.quit()
                # if Rsha_magma.round() == sha:
                #     print("-- agrees with Magma")
                # else:
                #     print("******** disagrees with Magma, which has {} =~= {}".format(Rsha_magma, Rsha_magma.round()))

            else:
                if verbose:
                    print("Unable to compute regulator or analytic Sha, since analytic rank = {} but we only have {} generators".format(ar, Edata['ngens']))
                Edata['sha'] = None

            line = make_mwdata_line(Edata)
            if verbose:
                print("New mwdata line: {}".format(line))
            mwdata.write(line + "\n")

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

