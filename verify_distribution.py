"""Check every file of the public source-only distribution (standard library).

Run python3 verify_distribution.py after unpacking the ancillary archive.
No historical manifests, exporter, author review or TeX files are required.
The check proves integrity and serialization, not the mathematical theorem.
"""
import argparse,gzip,hashlib,json,signal,time
from pathlib import Path
from read_coefficients import load
HERE=Path(__file__).resolve().parent
def main():
    start=time.monotonic()
    if hasattr(signal,'alarm'):signal.alarm(35)
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,default=HERE/'PUBLIC_MANIFEST.json')
    a=p.parse_args();manifest=json.loads(a.manifest.read_text());total=0
    for name,rec in manifest['files'].items():
        part=Path(name);assert not part.is_absolute() and '..' not in part.parts
        path=a.manifest.parent/part;assert path.resolve().is_relative_to(a.manifest.parent.resolve())
        b=path.read_bytes();assert len(b)==rec['bytes'] and hashlib.sha256(b).hexdigest()==rec['sha256'],name
        if name.endswith('.gz'):
            raw=gzip.decompress(b);assert len(raw)==rec['uncompressed_bytes'] and hashlib.sha256(raw).hexdigest()==rec['uncompressed_sha256']
            if name.endswith('.json.gz'):json.loads(raw)
        elif name.endswith('.json'):json.loads(b)
        total+=len(b)
    assert total==manifest['total_file_bytes'];z=load(a.manifest.parent/'f56_Z.txt.gz')
    print(json.dumps({'status':'PASS_PUBLIC_DISTRIBUTION_BYTES_AND_COEFFICIENT_ENCODING','files_checked':len(manifest['files']),'total_file_bytes':total,'polynomial_coefficients':len(z),'mathematical_proof_checked':False,'seconds':time.monotonic()-start},indent=2))
if __name__=='__main__':main()
