fields = {}
field_types = {2: 'RQF', -2:'IQF', 3: 'cubics', 4: 'quartics', 5: 'quintics', 6: 'sextics'}
dirs = {}

forms = None
nfcurves = None # keep pyflakes happy: the function below needs
                # nfcurves and forms to be set to the nfcurves collection.

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
    #ncu_CM = nfcurves.count({'field_label':fld, 'label': {'$regex':'CM'}})
    ncl_CM = nfcurves.find({'field_label':fld, 'number':int(1), 'label': {'$regex':'CM'}}).count()
    #ncu_nonCM = ncu-ncu_CM
    ncl_nonCM = ncl-ncl_CM
    if nforms!=ncl_nonCM:
        print("Field {} has {} rational newforms but {} non-CM isogeny classes".format(fld,nforms,ncl_nonCM))
        if nforms>ncl_nonCM:
            print("{} missing isogeny classes:".format(nforms-ncl_nonCM))
            for f in forms.find({'field_label':fld, 'dimension':int(1)}):
                if not nfcurves.find_one({'class_label':f['label']}):
                    print("Form {} has no matching curve".format(f['label']))
    else:
        print("Field {} has {} rational newforms and non-CM isogeny classes".format(fld,nforms))

    iso_file = "{}/isoclass.{}".format(pre,fld)
    cur_file = "{}/curves.{}".format(pre,fld)
    try:
        n_iso = len(open(iso_file).readlines())
    except:
        print("No file {} exists".format(iso_file))
        n_iso = 0

    if n_iso!= ncl:
        print("Field {}: file has {} classes, database has {}".format(fld,n_iso,ncl))
    else:
        print("Field {}: file and database both have {} classes".format(fld,ncl))

    try:
        n_cur = len(open(cur_file).readlines())
    except:
        print("No file %s exists" % cur_file)
        n_cur = 0

    if n_cur!= ncu:
        print("Field {}: file has {} curves, database has {}".format(fld,n_cur,ncu))
    else:
        print("Field {}: file and database both have {} curves".format(fld,ncu))

def check_data2(fld):
    """Wrapper round the find_curve_labels() function from
    hmf_check_find, which will check that for each rational newform
    there is a curve with the correct label, conductor and ap.
    """
    print("Field %s" % fld)
    find_curve_labels(fld)
