from lmfdb import db

DD=[15, 20, 24, 35, 39, 40, 51, 52, 55, 56, 68, 84, 87, 88, 91, 95]
dd=[15, 5, 6, 35, 39, 10, 51, 13, 55, 14, 17, 21, 87, 22, 91, 95]

def fix(D):
    # replace the curves:

    field_label = f"2.0.{D}.1"
    nc = db.ec_nfcurves.count({'field_label': field_label})
    print(f"{nc} curves exist for field {field_label}")
    print(" - deleting these...")
    db.ec_nfcurves.delete({'field_label': field_label}, restat=False)
    curve_file = f"/scratch/home/jcremona/ecnf-upload/ec_nfcurves.{field_label}"
    print(f" - deletion done, now reuploading from {curve_file}...")
    db.ec_nfcurves.copy_from(curve_file, restat=False, resort=False)
    ncnew = db.ec_nfcurves.count({'field_label': field_label})
    print(f" - done, now {ncnew} curves exist for field {field_label}")
    if nc!=ncnew:
       print(f" *********** numbers are not the same!")

    # replace the forms (no need to replace the dims):
    d = D if D%2 else D//4
    N1 = 1
    N2 = 1000
    forms_file = f"/scratch/home/jcremona/bmf-upload/bmf_forms.{d}.{N1}-{N2}"
    nf = db.bmf_forms.count({'field_label': field_label})
    print(f"{nf} forms exist for field {field_label}")
    print(" - deleting these...")
    db.bmf_forms.delete({'field_label': field_label}, restat=False)
    print(f" - deletion done, now reuploading from {forms_file}...")
    db.bmf_forms.copy_from(forms_file, restat=False, resort=False)
    nfnew = db.bmf_forms.count({'field_label': field_label})
    print(f" - done, now {nfnew} forms exist for field {field_label}")
    if nf!=nfnew:
       print(f" *********** numbers are not the same!")

