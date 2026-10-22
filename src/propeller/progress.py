"""Read-only inspection of durable runs, including an actively growing log."""
import json
from pathlib import Path
import re


def tail(path,limit=24000):
    if not path.exists():return ""
    with path.open('rb') as stream:
        stream.seek(0,2);size=stream.tell();stream.seek(max(0,size-limit))
        return stream.read().decode(errors='replace')


def run_status(root):
    records=[]
    for state_path in sorted(Path(root).glob('*/*/state.json'),key=lambda p:p.stat().st_mtime,reverse=True):
        state=json.loads(state_path.read_text());directory=state_path.parent
        record={"attempt":directory.name,"path":str(directory.resolve()),"status":state['status']}
        logs=list((directory/'case').glob('log.*'))
        if logs:
            log=max(logs,key=lambda p:p.stat().st_mtime);text=tail(log)
            record['latest_stage']=log.name.removeprefix('log.')
            iterations=re.findall(r'^Time = (\S+)',text,re.MULTILINE)
            if iterations:record['iteration']=iterations[-1]
            record['log_modified']=log.stat().st_mtime
        result=directory/'result.json'
        if result.exists():
            output=json.loads(result.read_text());record['reasons']=output.get('reasons',[])
            record['metrics']=output.get('metrics')
        records.append(record)
    return records
