fields = {}
field_types = {2: 'RQF', -2:'IQF', 3: 'cubics', 4: 'quartics', 5: 'quintics', 6: 'sextics'}
dirs = {}

def set_fields(degree=2):
    global fields, dirs
    if degree in fields:
        return
    fields[degree] = nfcurves.find({'degree':int(degree)}).distinct('field_label')
    dirs[degree] = "/home/jec/ecnf-data/{}".format(field_types[degree])


def check_data1(fld, pre, verbose=True):
    """Check that the number of rational newforms equals the number of
    elliptic curve isogeny classes in the database, and that the
    number of curves and classes in the database agree with the
    numbers of lines in the upload data files, for a given field
    label.

    ``pre`` is the directory where the relevant curves and isoclass
    data files are.

    """
    if verbose:
        print("Checking data for field {}".format(fld))
    nforms = forms.find({'field_label':fld, 'dimension':int(1)}).count()
    if verbose:
        print("{} rational newforms".format(nforms))
    if nforms==0:
        return
    ncu = nfcurves.find({'field_label':fld}).count()
    ncl = nfcurves.find({'field_label':fld, 'number':int(1)}).count()
    if nforms!=ncl:
        print("Field %s has %s rational newforms but %s isogeny classes" % (fld,nforms,ncl))

    iso_file = "%s/isoclass.%s" % (pre,fld)
    cur_file = "%s/curves.%s" % (pre,fld)
    try:
        n_iso = len(file(iso_file).readlines())
    except:
        print("No file %s exists" % iso_file)
        n_iso = 0
    if n_iso!= ncl:
        print("Field %s: file has %s classes, database has %s" % (fld,n_iso,ncl))

    try:
        n_cur = len(file(cur_file).readlines())
    except:
        print("No file %s exists" % cur_file)
        n_cur = 0

    if n_cur!= ncu:
        print("Field %s: file has %s curves, database has %s" % (fld,n_cur,ncu))

def check_data2(fld):
    """Wrapper round the find_curve_labels() function from
    hmf_check_find, which will check that for each rational newform
    there is a curve with the correct label, conductor and ap.
    """
    print("Field %s" % fld)
    find_curve_labels(fld)
