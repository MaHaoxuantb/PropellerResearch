"""Content-addressed experiments with immutable attempts and durable results."""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sqlite3
import shutil
import time
import uuid


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False)


def source_fingerprint():
    root=Path(__file__).parent
    digest=hashlib.sha256()
    for path in sorted(root.rglob("*.py")):
        digest.update(str(path.relative_to(root)).encode());digest.update(path.read_bytes())
    return digest.hexdigest()


def versions():
    result={}
    for package in ("numpy","trimesh","manifold3d","scipy","networkx","shapely","rtree"):
        try: result[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: result[package]=None
    return result


def evaluation_key(config):
    # Include all inputs conservatively; never merge nearby float vectors.
    payload={"config":config,"source":source_fingerprint(),"dependencies":versions()}
    return hashlib.sha256(canonical(payload).encode()).hexdigest(),payload


def atomic_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+"."+uuid.uuid4().hex+".tmp")
    with temporary.open("w") as stream:
        stream.write(json.dumps(value,indent=2,allow_nan=False));stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)


class Store:
    def __init__(self,root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.connection=sqlite3.connect(self.root/"index.sqlite",timeout=30)
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("CREATE TABLE IF NOT EXISTS attempts (id TEXT PRIMARY KEY, cache_key TEXT, status TEXT, path TEXT, started REAL, finished REAL)")
        self.connection.commit()

    def recover_interrupted(self):
        """Mark only locally recorded attempts whose owner process is absent.

        A permission error is not evidence that a process is dead. Existing
        immutable results are never rewritten by recovery.
        """
        recovered=[]
        for identity,directory in self.connection.execute("SELECT id,path FROM attempts WHERE status='RUNNING'").fetchall():
            state_path=Path(directory)/"state.json"
            if not state_path.exists():continue
            state=json.loads(state_path.read_text());pid=state.get("pid")
            if not isinstance(pid,int) or pid<=0:continue
            try:os.kill(pid,0)
            except ProcessLookupError:
                if not (Path(directory)/"result.json").exists():
                    self.finish(identity,Path(directory),{"status":"INTERRUPTED","reasons":["Evaluation owner process exited without publishing a result"]})
                    recovered.append(identity)
            except PermissionError:continue
        return recovered

    def cached(self,key):
        row=self.connection.execute("SELECT path FROM attempts WHERE cache_key=? AND status='SUCCESS' ORDER BY finished DESC LIMIT 1",(key,)).fetchone()
        if row:
            file=Path(row[0])/"result.json"
            if file.exists(): return json.loads(file.read_text())
        return None

    def start(self,key,manifest):
        identity=uuid.uuid4().hex
        directory=(self.root/key[:16]/identity).resolve();directory.mkdir(parents=True)
        atomic_json(directory/"manifest.json",manifest)
        source=Path(__file__).parent
        for path in source.rglob("*.py"):
            destination=directory/"source"/"propeller"/path.relative_to(source)
            destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(path,destination)
        atomic_json(directory/"state.json",{"status":"RUNNING","pid":os.getpid(),"started":time.time()})
        self.connection.execute("INSERT INTO attempts VALUES (?,?,?,?,?,?)",(identity,key,"RUNNING",str(directory),time.time(),None));self.connection.commit()
        return identity,directory

    def finish(self,identity,directory,result):
        result=dict(result,attempt_id=identity,path=str(directory))
        atomic_json(directory/"result.json",result)
        atomic_json(directory/"state.json",{"status":result["status"],"finished":time.time()})
        self.connection.execute("UPDATE attempts SET status=?, finished=? WHERE id=?",(result["status"],time.time(),identity));self.connection.commit()
        return result

    def close(self): self.connection.close()
