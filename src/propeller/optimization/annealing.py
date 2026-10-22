"""Bounded SA with separate feasible archive and exact RNG checkpoint resume."""
import copy
import math
from pathlib import Path
import json
import time
import numpy as np
from ..config import normalized,physical,with_parameters
from ..experiments import canonical,atomic_json,source_fingerprint


def reflect(q):
    q=np.mod(q,2)
    return np.where(q>1,2-q,q)


def score(result,config):
    if result.get("status")!="SUCCESS": return None
    metrics=result.get("metrics",{});T,P=metrics.get("thrust"),metrics.get("shaft_power")
    if T is None or P is None or not math.isfinite(T) or not math.isfinite(P) or P<=0: return None
    o=config["operating"];target=o["target_thrust"];reference=o["reference_power"]
    F=P/reference+config["optimization"]["thrust_penalty"]*max(0,1-T/target)**2
    return {"score":F,"feasible":T>=target,"thrust":T,"power":P}


def optimize(config,evaluator,checkpoint,resume=False,proposal_limit=None):
    o=config["operating"];settings=config["optimization"]
    if o["target_thrust"] is None or o["reference_power"] is None:
        raise ValueError("Freeze positive target_thrust and reference_power from a converged baseline before optimization")
    checkpoint=Path(checkpoint);rng=np.random.default_rng(settings["seed"])
    fingerprint=canonical({"config":config,"source":source_fingerprint()})
    if resume:
        state=json.loads(checkpoint.read_text())
        if state["fingerprint"]!=fingerprint: raise ValueError("Checkpoint configuration/source differs")
        rng.bit_generator.state=state["rng"]
    else:
        if checkpoint.exists(): raise FileExistsError("Checkpoint exists; use resume or a new path")
        result=evaluator(config);current=score(result,config)
        if current is None: raise ValueError("Initial candidate has no successful finite aerodynamic score")
        state={"fingerprint":fingerprint,"q":normalized(config),"current":current,"best":None,
               "temperature":settings["temperature"],"proposals":0,"evaluations":0 if result.get("cache_hit") else 1,
               "scored":1,"history":[],"elapsed":0.0}
        if current["feasible"]: state["best"]={"parameters":config["parameters"],**current}
    start=time.monotonic();initial_proposals=state["proposals"]
    def save():
        state["rng"]=rng.bit_generator.state
        atomic_json(checkpoint,state)
    save()
    try:
        while state["proposals"]<settings["max_proposals"] and state["evaluations"]<settings["max_evaluations"]:
            if state["elapsed"]+time.monotonic()-start>=settings["max_seconds"]: break
            if proposal_limit is not None and state["proposals"]-initial_proposals>=proposal_limit: break
            # Persist proposed point before expensive evaluation, so a restart replays it.
            if "pending" in state:
                q=np.array(state["pending"])
            else:
                q=np.array(state["q"]);step=np.zeros(len(q))
                if rng.random()<.75:
                    index=rng.integers(len(q));step[index]=rng.normal(0,settings["proposal_scale"])
                else: step=rng.normal(0,settings["proposal_scale"],len(q))
                q=reflect(q+step);state["pending"]=q.tolist();save()
            try:
                candidate=with_parameters(config,physical(config,q))
                result=evaluator(candidate)
            except ValueError as exc: result={"status":"GEOMETRY_INVALID","reasons":[str(exc)],"cache_hit":True}
            value=score(result,config);accepted=False
            if not result.get("cache_hit",False) and result.get("status")!="GEOMETRY_INVALID":state["evaluations"]+=1
            if value is not None:
                delta=value["score"]-state["current"]["score"]
                accepted=delta<=0 or rng.random()<math.exp(-delta/max(state["temperature"],1e-15))
                if accepted:state["q"]=q.tolist();state["current"]=value
                if value["feasible"] and (state["best"] is None or value["power"]<state["best"]["power"]):
                    state["best"]={"parameters":physical(config,q),**value}
                state["temperature"]*=settings["cooling"];state["scored"]+=1
            state["proposals"]+=1
            state["history"].append({"proposal":state["proposals"],"parameters":physical(config,q),"status":result["status"],"accepted":accepted,"value":value,"attempt_id":result.get("attempt_id"),"cache_hit":result.get("cache_hit",False)})
            del state["pending"];save()
    finally:
        state["elapsed"]+=time.monotonic()-start;save()
    return state
