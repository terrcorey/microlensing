"""Download raw photometry for each event we're working with.

Sources are redistributed (MIT license) by the MulensModel project for
tutorial use, originally from OGLE/MOA. OGLE's own EWS terms ask to be
contacted/credited before publishing anything built on their data.
"""

import urllib.request
from pathlib import Path

BASE = "https://raw.githubusercontent.com/rpoleski/MulensModel/master/data/photometry_files"

SOURCES = {
    # single point-lens event, see dataset_names.txt (O-05-BLG086)
    "OGLE-2005-BLG-086.dat": f"{BASE}/OB05086/starBLG234.6.I.218982.dat",
    # binary-lens planet-discovery event (Bond et al. 2004), see dataset_names.txt (O-03-BLG235)
    "OGLE-2003-BLG-235_OGLE.tbl.txt": f"{BASE}/OB03235/OB03235_OGLE.tbl.txt",
    "OGLE-2003-BLG-235_MOA.tbl.txt": f"{BASE}/OB03235/OB03235_MOA.tbl.txt",
}

DATA_DIR = Path(__file__).parent / "data"

if __name__ == "__main__":
    DATA_DIR.mkdir(exist_ok=True)
    for filename, url in SOURCES.items():
        out_path = DATA_DIR / filename
        if out_path.exists():
            print(f"already have {out_path}")
            continue
        urllib.request.urlretrieve(url, out_path)
        print(f"downloaded {out_path}")
