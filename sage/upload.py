import os
import sys
from sage.all import Set, ZZ, RR
from files import ECNF_DIR, ECNF_UPLOAD_DIR, all_ftypes, read_all_field_data
from schemas import ec_nfcurves_schema, keys_and_types, ec_nfcurves_extra_columns, extra_keys_and_types

# Any columns in ec_nfcurves which are lists must be in this list so
# that the upload file uses "{...}" instead of "[...]":

postgres_array_cols = [k for k,v in ec_nfcurves_schema.items() if '[' in v]

# ['nonmax_primes', 'heights', 'isodeg', 'torsion_primes', 'reducible_primes', 'conductor_norm_factors']


def get_column_names_and_types():
    t = ec_nfcurves_schema
    return {k: t[k] for k in sorted(t)}

ec_nfcurves_column_names_and_types = get_column_names_and_types()
ec_nfcurves_columns = Set(ec_nfcurves_column_names_and_types.keys())

assert ec_nfcurves_columns == Set(keys_and_types.keys()) + Set(['id'])
assert Set(ec_nfcurves_extra_columns)==Set(extra_keys_and_types.keys())

def column_to_string(colname, col):
        if col is None:
            return "\\N"
        col = str(col).replace(" ","")
        col = col.replace("'",'"')
        if col == "True":
            return "t"
        if col == "False":
            return "f"
        if colname in postgres_array_cols:
            col = col.replace("[","{").replace("]","}")
        if colname == 'equation':
            col = col.replace("\\","\\\\")
            #col = ''.join(['"', col, '"'])
        if colname == 'local_data':
            col = col.replace("None", "null")
        if colname == 'disc' and "." in col:
            print("Old disc: {}".format(col))
            col = "({})".format(ZZ(RR(col[1:-1])))
            print("New disc: {}".format(col))
        return col

def data_to_string(n, record, columns=None):
    """NB A list stored in the database as a postgres array (e.g. int[]
    or numeric[]) must appear as (e.g.) {1,2,3} not [1,2,3].

    The relevant columns are in the global array postgres_array_cols

    If columns is not None, then only these columns are output.
    """
    record['id'] = n
    if columns:
        keys = columns
        if 'label' in keys:
            keys.remove('label')
    else:
        keys = list(ec_nfcurves_column_names_and_types.keys())
        keys.remove('id')
        keys.remove('label')
    keys = ['label'] + keys
    return "|".join([column_to_string(k, record[k]) for k in keys])

def make_upload_file(ftypes=all_ftypes, fields=None, xfields=None, columns=None, outfilename=None):
    """
    ftypes: list of one or more from 'IQF', 'RQF', 'cubics', 'gunnells', 'quartics', 'quintics', 'sextics'
            (default: all of these field types are processed)

    fields: list of one or more field labels; only makes sense if ftypes has one entry.
            (default: all fields of the field types requested are processed)

    xfields: list of one or more field labels to omit; only makes sense if ftypes has one entry.
            (default: all fields of the field types requested are processed)

    columns: if None (default), output files contains all columns.  Otherwise just output these columns.

    outfilesname: if None (default) output to stdout, else output to this file (which will be overwritten if it exists)
    """
    alldata = {}
    for ftype in ftypes:
        print("reading data from {}".format(ftype))
        if not fields:
            with open("{}/{}/fields.txt".format(ECNF_DIR,ftype)) as field_file:
                flds = [f[:-1] for f in field_file.readlines()]
        else:
            flds = fields
        if xfields:
            for f in xfields:
                print("excluding field {}".format(f))
                flds.remove(f)
        for fname in flds:
            print("reading data for {}".format(fname))
            data = read_all_field_data(os.path.join(ECNF_DIR,ftype), fname, check_cols=True, mwdata_format="new")
            alldata.update(data)
    if outfilename:
        outfile = open(os.path.join(ECNF_UPLOAD_DIR,outfilename), 'w')
    else:
        outfile = sys.stdout
    if columns:
        keys = columns
        if 'label' in keys:
            keys.remove('label')
    else:
        keys = list(ec_nfcurves_column_names_and_types.keys())
        keys.remove('id')
        keys.remove('label')
    keys = ['label'] + keys
    vals = [ec_nfcurves_column_names_and_types[k] for k in keys]
    outfile.write("|".join(keys))
    outfile.write("\n")
    outfile.write("|".join(vals))
    outfile.write("\n\n")
    id = 0
    for label in sorted(alldata):
        id += 1
        alldata[label]['non-surjective_primes'] = alldata[label]['nonmax_primes']
        outfile.write(data_to_string(id, alldata[label], columns=columns))
        outfile.write("\n")
    if outfilename:
        outfile.close()
    print("{} data lines + 3 header lines with {} columns written to {}".format(id,len(keys),outfilename if outfilename else "stdout"))
    print("Columns:\n{}".format(keys))
