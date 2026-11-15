"""Verify bibliographic records against the DOI registration agency."""
import concurrent.futures
import json
from pathlib import Path
import urllib.request
import time

DOIS={
 'stahl2016':'10.1126/science.aaf2403',
 'rodriques2019':'10.1126/science.aaw1219',
 'moncada2020':'10.1038/s41587-019-0392-8',
 'dries2021':'10.1186/s13059-021-02286-2',
 'palla2022':'10.1038/s41592-021-01358-2',
 'hao2021':'10.1016/j.cell.2021.04.048',
 'wolf2018':'10.1186/s13059-017-1382-0',
 'cable2022':'10.1038/s41587-021-00830-w',
 'kleshchevnikov2022':'10.1038/s41587-021-01139-4',
 'elosua2021':'10.1093/nar/gkab043',
 'hanzelmann2013':'10.1186/1471-2105-14-7',
 'andreatta2021':'10.1016/j.csbj.2021.06.043',
 'phipson2010':'10.2202/1544-6115.1585',
 'bh1995':'10.1111/j.2517-6161.1995.tb02031.x',
 'virtanen2020':'10.1038/s41592-019-0686-2',
 'sandve2013':'10.1371/journal.pcbi.1003285',
 'wilkinson2016':'10.1038/sdata.2016.18',
 'liu2025':'10.1016/j.ccell.2025.03.004',
 'alam2025':'10.1101/2025.05.08.652741',
 'zhong2024':'10.1038/s41592-024-02463-8',
 'tang2026':'10.1038/s41556-025-01838-z',
 'lahiri2006':'10.1214/009053606000000551',
}


def fetch(item):
    key,doi=item
    req=urllib.request.Request('https://api.crossref.org/works/'+doi,headers={'User-Agent':'CAD-SPOT bibliography verification/2.0'})
    try:
        with urllib.request.urlopen(req,timeout=40) as response:
            data=json.load(response)['message']
        return key,dict(doi=doi,status='verified_crossref',record=data)
    except Exception as e:
        return key,dict(doi=doi,status='unverified',error=str(e))


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    path=root/'Review/verified_references.json'
    data=json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    for key,doi in DOIS.items():
        if data.get(key,{}).get('status')=='verified_crossref':
            continue
        for attempt in range(3):
            _,value=fetch((key,doi))
            data[key]=value
            if value['status']=='verified_crossref':
                break
            time.sleep(2+attempt)
        time.sleep(1)
    (root/'Review/verified_references.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for key,value in data.items():
        print(key,value['status'],value.get('record',{}).get('title',[''])[0],flush=True)
