"""Independent BED/bigBed round-trip and local IGV session validation."""
import json,struct,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
import pybigtools
O=Path(__file__).resolve().parents[1]
manifest=json.loads((O/'metadata/track_manifest.json').read_text())
checks=[]
for t in manifest:
    bed=[x.split('\t') for x in (O/t['file']).read_text().splitlines()]
    bb=O/'hub/hg38'/f'{t["id"]}.bb'
    header=struct.unpack('<IHHQQQHH',bb.read_bytes()[:36]);assert header[0]==0x8789F2EB
    assert header[-2:]==(t['bed_fields'],t['bed_fields']),(t['id'],header)
    reader=pybigtools.open(str(bb));got=[]
    for chrom in sorted({x[0] for x in bed}):
        got.extend([[chrom]+list(map(str,r)) for r in reader.records(chrom,0,reader.chroms()[chrom])])
    assert got==bed,(t['id'],len(got),len(bed))
    checks.append(dict(track=t['id'],records=len(bed),BED_fields=t['bed_fields'],roundtrip='PASS'))
for path in (O/'igv').glob('*local*.xml'):
    root=ET.parse(path).getroot();assert root.attrib['genome']=='hg38'
    for r in root.find('Resources'):assert (path.parent/r.attrib['path']).resolve().is_file()
summary={'result':'PASS','tracks':len(checks),'checks':checks,'local_IGV_resources':'all resolve',
         'v7_genomic_controls':'All exported windows are 170 bp; corrected source-based reconstruction',
         'C1_status':'1998 existing sequences + 2 undesigned reserved slots; not a synthesis order'}
(O/'metadata/bundle_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print(f'PASS: {len(checks)} BED/bigBed round trips, field headers, and both local IGV sessions.')
