"""Download raw photometry for each event we're working with.

O-05-BLG086/O-03-BLG235 are redistributed (MIT license) by the MulensModel
project for tutorial use, originally from OGLE/MOA; O-05-BLG169 comes from the
NASA Exoplanet Archive. OGLE's own EWS terms ask to be contacted/credited
before publishing anything built on their data.
"""

import urllib.request
from pathlib import Path

MM_BASE = "https://raw.githubusercontent.com/rpoleski/MulensModel/master/data/photometry_files"
NASA_BASE = "https://exoplanetarchive.ipac.caltech.edu/data/ExoData/0300/0300030/data"

SOURCES = {
    # single point-lens event, see dataset_names.txt (O-05-BLG086)
    "OGLE-2005-BLG-086.dat": f"{MM_BASE}/OB05086/starBLG234.6.I.218982.dat",
    # binary-lens planet-discovery event (Bond et al. 2004), see dataset_names.txt (O-03-BLG235)
    "OGLE-2003-BLG-235_OGLE.tbl.txt": f"{MM_BASE}/OB03235/OB03235_OGLE.tbl.txt",
    "OGLE-2003-BLG-235_MOA.tbl.txt": f"{MM_BASE}/OB03235/OB03235_MOA.tbl.txt",
    # binary-lens planet event (Gould et al. 2006, archive UID 0300030), see dataset_names.txt (O-05-BLG169)
    "OGLE-2005-BLG-169_Auckland.tbl": f"{NASA_BASE}/UID_0300030_PLC_001.tbl",  # unfiltered
    "OGLE-2005-BLG-169_FTN.tbl": f"{NASA_BASE}/UID_0300030_PLC_002.tbl",  # R
    "OGLE-2005-BLG-169_MDM.tbl": f"{NASA_BASE}/UID_0300030_PLC_003.tbl",  # I
    "OGLE-2005-BLG-169_OGLE.tbl": f"{NASA_BASE}/UID_0300030_PLC_004.tbl",  # I
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
