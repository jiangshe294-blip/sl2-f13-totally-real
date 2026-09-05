"""Read the actual degree-56 polynomial without million-digit Python integers.

Python 3.8+, standard library only. A rational printed as n/D is exact;
no cancellation is necessary to specify that rational number.
"""
import argparse, gzip, hashlib, json, re
from pathlib import Path

HERE=Path(__file__).resolve().parent
PACKED_SHA='82343cb0db291f1b519a418957ac00f455d9dc07d27d58c84d0245517950d113'
TEXT_SHA='5058792bf79dd594034393954aac53bf801a57f7a36e989c55efefc9c270fd50'
def load(path=None):
    path=Path(path) if path else HERE/'f56_Z.txt.gz'
    packed=path.read_bytes()
    assert hashlib.sha256(packed).hexdigest()==PACKED_SHA,'gzip digest mismatch'
    text=gzip.decompress(packed)
    assert hashlib.sha256(text).hexdigest()==TEXT_SHA,'coefficient-text digest mismatch'
    assert text.endswith(b'\n')
    values=text.decode('ascii').splitlines()
    assert len(values)==57
    assert all(re.fullmatch(r'0|-?[1-9][0-9]*',x) for x in values)
    assert all(values[i]=='0' for i in range(1,57,2))
    assert sum(x!='0' for x in values)==29
    assert values[-1][0] in '123456789'
    return values
def rational(z,D):
    if z=='0':return '0'
    if z==D:return '1'
    return z+'/'+D
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--path',type=Path)
    p.add_argument('--form',choices=['Z','Q','P28'],default='Z')
    p.add_argument('--coefficient',type=int)
    p.add_argument('--all',action='store_true')
    a=p.parse_args();z=load(a.path);D=z[-1]
    values=z if a.form=='Z' else [rational(v,D) for v in z]
    if a.form=='P28':values=values[::2]
    if a.coefficient is not None:
        if not 0<=a.coefficient<len(values):p.error('coefficient outside polynomial degree')
        print(values[a.coefficient])
    elif a.all:print('\n'.join(values))
    else:print(json.dumps({'degree':56,'ascending_coefficients':57,'nonzero_coefficients':29,'D_decimal_digits':len(D),'f_definition':'sum z_i X^i / D, D=z_56','P28_definition':'sum z_(2i) T^i / D','compressed_sha256':PACKED_SHA,'uncompressed_sha256':TEXT_SHA,'checks':'hashes, canonical integer spelling, count, parity and positive leading coefficient; content and algebraic proof checked separately'},indent=2))
if __name__=='__main__':main()
