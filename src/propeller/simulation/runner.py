"""Native launcher adapter. Commands are argument lists, never shell input."""
import os
import hashlib
from pathlib import Path
import signal
import subprocess
import time
import resource
import sys
from ..experiments import atomic_json


def process_alive(pid):
    if not isinstance(pid,int) or pid<=0: return False
    try: os.kill(pid,0)
    except ProcessLookupError: return False
    except PermissionError: return True
    return True


def runtime_info(config):
    launcher=Path(config["cfd"]["launcher"])
    if not launcher.is_file(): return {"available":False,"reason":"Native launcher not found","launcher":str(launcher)}
    try:
        # Trusted constant script; no config values interpolated into shell text.
        r=subprocess.run([str(launcher),"bash","-c",'printf "%s\\n" "$WM_PROJECT_VERSION"; command -v simpleFoam blockMesh snappyHexMesh'],capture_output=True,text=True,timeout=30)
        lines=r.stdout.strip().splitlines()
        binaries={}
        for entry in lines[1:]:
            file=Path(entry)
            if file.is_file():binaries[str(file)]=hashlib.sha256(file.read_bytes()).hexdigest()
        if len(lines)>1:
            library=Path(lines[1]).parent.parent/"lib"
            for name in ("libOpenFOAM.dylib","libfiniteVolume.dylib"):
                file=library/name
                if file.is_file():binaries[str(file)]=hashlib.sha256(file.read_bytes()).hexdigest()
        return {"available":r.returncode==0 and len(lines)>=4 and lines[0]==config["cfd"]["version"],
                "version":lines[0] if lines else None,"executables":lines[1:],"binary_sha256":binaries,"launcher":str(launcher),"stderr":r.stderr}
    except (OSError,subprocess.TimeoutExpired) as exc: return {"available":False,"reason":str(exc)}


def run_command(config,case,command):
    label="simpleFoam" if command[0]=="mpirun" else command[0]
    logfile=Path(case)/("log."+label)
    args=[config["cfd"]["launcher"]]+list(command)+["-case",str(Path(case).resolve())]
    started=time.monotonic()
    absolute_deadline=float(os.environ.get('PROPELLER_DEADLINE_EPOCH','inf'))
    supervisor=int(os.environ.get('PROPELLER_SUPERVISOR_PID','0'))
    if time.time()>=absolute_deadline or (supervisor and not process_alive(supervisor)):
        return {'command':args,'returncode':None,'status':'TIMEOUT','seconds':0,'log':str(logfile),'reason':'campaign deadline or supervisor absent'}
    with logfile.open("w") as log:
        process=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        atomic_json(Path(case)/'active-process.json',{'pid':process.pid,'owner_pid':os.getpid(),'command':args,'started':time.time(),'status':'RUNNING'})
        try:
            deadline=started+config["cfd"]["timeout_seconds"]
            startup=started+config["cfd"].get("startup_timeout_seconds",60)
            while process.poll() is None:
                if time.monotonic()>=deadline or time.time()>=absolute_deadline or (supervisor and not process_alive(supervisor)) or (time.monotonic()>=startup and logfile.stat().st_size==0):
                    raise subprocess.TimeoutExpired(args,config["cfd"]["timeout_seconds"])
                time.sleep(.2)
            code=process.returncode
        except (subprocess.TimeoutExpired,KeyboardInterrupt) as exc:
            try:os.killpg(process.pid,signal.SIGTERM)
            except ProcessLookupError:pass
            try: process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait()
            if isinstance(exc,KeyboardInterrupt):
                atomic_json(Path(case)/'active-process.json',{'pid':process.pid,'status':'INTERRUPTED'})
                raise
            atomic_json(Path(case)/'active-process.json',{'pid':process.pid,'status':'TIMEOUT'})
            return {"command":args,"returncode":None,"status":"TIMEOUT","seconds":time.monotonic()-started,"log":str(logfile)}
    atomic_json(Path(case)/'active-process.json',{'pid':process.pid,'status':'FINISHED','returncode':code})
    rss=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return {"command":args,"returncode":code,"status":"SUCCESS" if code==0 else "FAILED","seconds":time.monotonic()-started,"log":str(logfile),
            "cumulative_child_max_rss_bytes":int(rss if sys.platform=='darwin' else rss*1024)}
