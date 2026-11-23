"""Read local OPS description XML without flattening paragraph/table provenance."""
from pathlib import Path
import html
import xml.etree.ElementTree as ET

def read_description(path):
    path=Path(path);root=ET.parse(path).getroot();segments=[]
    local=lambda e:e.tag.rsplit('}',1)[-1]
    descriptions=[e for e in root.iter() if local(e)=='description']
    if not descriptions:raise ValueError('No description element in EPO XML')
    for description in descriptions:
        for e in description:
            raw=ET.tostring(e,encoding='unicode');text=html.unescape(''.join(e.itertext()))
            segments.append({'index':len(segments),'tag':local(e),'id':e.get('id'),'num':e.get('num'),
                             'text':text,'xml':raw,'has_image':('<img' in text.lower() or any(local(x) in ('img','image') for x in e.iter()))})
    return {'source':str(path.resolve()),'source_type':'epo_description','indexes_are_pdf_pages':False,'segments':segments}
