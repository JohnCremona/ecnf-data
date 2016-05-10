quadratics = [f['label'] for f in fields.find({'degree':int(2)})]
print("%s quadratic fields" % len(quadratics))
pre = "/home/jec/ecnf-data/RQF"

# cubics     = [f['label'] for f in fields.find({'degree':int(3)})]
# print("%s cubic fields" % len(cubics))
# pre = "/home/jec/ecnf-data/cubics"

# quartics   = [f['label'] for f in fields.find({'degree':int(4)})]
# print("%s quartic fields" % len(quartics))
# pre = "/home/jec/ecnf-data/quartics"

def check_data1(fld, pre):
    """Check that the number of rational newforms equals the number if
    elliptic curve isogeny classes in the database, and that the
    number of curves and classes in the database agree with the
    numbers of lines in the upload data files, for a given field
    label.

    ``pre`` is the directory where the relevant curves and isoclass
    data files are.

    """
    nforms = forms.find({'field_label':fld, 'dimension':int(1)}).count()
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
