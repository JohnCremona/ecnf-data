# Sage interface to Sutherland's Magma script for Galois images

import re
import os
from sage.all import prod

HOME = os.getenv("HOME")
GALREP_SCRIPT_DIR = os.path.join(HOME, "galrep")

def init_galrep(mag, script_dir=GALREP_SCRIPT_DIR):
    """
    Load the 2adic magma script into this magma process
    """
    mag.eval('cwd:=GetCurrentDirectory();')
    mag.eval('ChangeDirectory("{}");'.format(script_dir))
    mag.eval('load "nfgalrep.m";')
    mag.eval('ChangeDirectory(cwd);')

def split_galois_image_code(s):
    """Each code starts with a prime (1-3 digits but we allow for more)
    followed by an image code for that prime.  This function returns
    two substrings, the prefix number and the rest.
    """
    p = re.findall(r'\d+', s)[0]
    return p, s[len(p):]

def galrep_data_from_magma(E, mag):
    """Use Magma script to compute mod-p Galois image data

    E is an elliptic curve over a number field

    mag is a magma process.

    NB before calling this, the caller must have called init_galrep()

    Returns a single space-separated string containing all image codes
    at non-maximal primes (possibly empty)

    """
    return str(mag.ComputeGaloisImage(E))

# NB for elliptic curves over Q and over number fields, the columns
# holding the mod-p Galois image data are different!
#
# EC/Q           EC/NF
#
# galois_images  galois_images  (the same: one, possibly empty, string)
# modp_images    -              (possibly empty list of strings)
# nonmax_primes  nonmax_primes  (list of non-maximal primes)
# nonmax_rad     nonmax_rad     (product of non-maximal primes)


def parse_galrep_data_string(galois_images, verbose=False):
    if verbose:
        print("galois images = {}".format(galois_images))
    image_codes = galois_images.split() # list of strings
    pr = [int(split_galois_image_code(s)[0]) for s in image_codes]
    record = {'galois_images': image_codes,
              #'modp_images': image_codes,
              'nonmax_primes': pr,
              'nonmax_rad': prod(pr),
             }
    if verbose:
        print("galrep data: {}".format(record))
    return record

def get_galrep_data(E, mag=None, verbose=False):
    """
    Use Magma script to compute mod-p Galois image data

    E is an elliptic curve over a number field

    mag is a magma process.

    NB before calling this, the caller must have called init_galrep()

    Returns a dict with keys

    'galois_images': one string, possible empty
    'modp_images': list of strings, one per non-maximal prime
    'nonmax_primes': list of non-maximal primes
    'nonmax_rad': product of non-maximal primes
    """
    if not mag:
        from magma import get_magma
        mag = get_magma()
    return parse_galrep_data_string(galrep_data_from_magma(E, mag), verbose=verbose)


