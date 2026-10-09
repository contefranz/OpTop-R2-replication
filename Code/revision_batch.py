#!/usr/bin/env python3
"""Evaluation-only revision batch (the "last passage").

No topic model is estimated: every production stage runs with OPTOP_NO_FIT=1, so
a fit-cache miss is an error, and the sha256 of every file in Data/FITS is
compared with the manifest recorded before the revision, before and after the
run. Nothing is overwritten: every output carries the suffix (default _rev2).

A stage is COMPLETE only if it exited 0, all its declared outputs exist, and its
identity -- command, input hashes, package, and the code that determines its
result -- is the one on disk now. Identity is per stage: the shared modules, the
configuration and the script(s) the stage runs; repairing a figure script, or a
different driver, therefore never invalidates a night of scoring.
Inside a stage, the R drivers keep their own checkpoints (by seed / battery /
evaluation block, Data/Checkpoints/<suffix>/), reused only on identical
configuration, package and scoring code.

    python3 Code/revision_batch.py --select A          # night 1
    python3 Code/revision_batch.py --select B          # nights 2-3
    python3 Code/revision_batch.py --select final      # post-processing + acceptance
    python3 Code/revision_batch.py --select all --dry-run

Stage ids 2-8 are those of the revision plan and of the % PENDING comments in the
manuscript; 9-11 are the additional analyses; 12-13 close the batch.
"""
import argparse
import csv
import hashlib
import json
import os
import signal
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)

BASE = "Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01"
PRIOR = BASE.replace("fa0p1", "fa0p5")
DGP2 = "Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01"
PREP = "Data/MDNA/mdna_prep_2015_2016.qs2"
MDNA_TESTS_K = "40,50,60,170,180,190"      # the verified diagnostic set
BLOCKS = {"A": [2, 3, 4, 5, 6, 9, 10, 11], "B": [7, 8], "final": [12, 13]}
ORDER = [0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 7, 8, 12, 13]   # execution order


def cache(exp, tag=BASE):
    return f"Data/{exp}/{exp.lower()}_results_{exp}_full_{tag}.qs2"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log(msg):
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {msg}", flush=True)


def save_json(path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def build_stages(w, sfx):
    """The stage table. `code`: files beyond the scoring core that determine the
    stage's result. The number of workers is NOT part of a stage's identity."""
    S = []

    def add(n, key, title, cmd, outputs=(), inputs=(), deps=(), code=(), eta="",
            blas="1"):
        # blas: BLAS/OpenMP threads of the stage's R process. Simulation stages
        # parallelise over forked workers, so each worker gets ONE thread (no
        # oversubscription). The MD&A stages score sequentially in the main
        # process on 3,466 x 27,576 dense products: they keep the library default,
        # which is also the threading the published MD&A numbers were computed with.
        S.append(dict(n=n, key=key, title=title, cmd=list(cmd), outputs=list(outputs),
                      inputs=list(inputs), deps=list(deps), code=list(code), eta=eta,
                      blas=blas))

    def rev(p):
        return p[:-4] + sfx + ".qs2"

    def sim(driver, src, *extra):
        return ["Rscript", f"Code/run_{driver}.R", "full", str(w),
                f"cfg_from={src}", f"out_suffix={sfx}", *extra]

    add(0, "units", "unit gates", ["Rscript", "Code/tests_unit.R"],
        code=["Code/tests_unit.R"], eta="<1 min")
    add(1, "smoke", "every modified driver end to end on toy corpora",
        ["bash", "Code/revision_smoke.sh", str(min(w, 2)), sfx],
        code=["Code/revision_smoke.sh", *sorted(str(p) for p in Path("Code").glob("run_*.R"))],
        eta="~3 min")

    def mdna(n, key, title, tag, c, grid, tests, eta):
        add(n, key, title,
            ["Rscript", "Code/run_mdna.R", f"workers={w}", f"K_grid={grid}",
             "refine_span=0", f"tests_K={tests}", f"c_part={c}", "min_null=1",
             f"out_suffix={tag}"],
            [f"Data/MDNA/mdna_results_MDNA_2015_2016{tag}.qs2"], [PREP], eta=eta,
            blas="default")

    mdna(2, "mdna", "MD&A full rescoring: word curves over K, raw word discrepancies, "
         "tests at the verified K set, support resolution", sfx, 1, "10:200:10",
         MDNA_TESTS_K, "~1.5 h")
    mdna(3, "mdna_c05", "MD&A design: c = 0.5, delta = 1, same 20-model grid",
         "_c05" + sfx, 0.5, "10:200:10", MDNA_TESTS_K, "~1.5 h")
    mdna(3, "mdna_c2", "MD&A design: c = 2, delta = 1, same 20-model grid",
         "_c2" + sfx, 2, "10:200:10", MDNA_TESTS_K, "~1.5 h")
    mdna(3, "mdna_g100", "MD&A design: grid truncated at K = 100 (c = 1, delta = 1)",
         "_g100" + sfx, 1, "10:100:10", "40,50,60", "~40 min")

    for name, tag in (("base", BASE), ("dgp2", DGP2)):
        src = cache("E2", tag)
        add(4, f"e2_{name}", f"E2 ({name}): 20,000 reference documents per training fit, "
            "reference SEs, three selection rules", sim("E2", src, "J_truth=20000"),
            [rev(src)], [src], eta="~1.5 h")
    for name, tag in (("base", BASE), ("dgp2", DGP2)):
        src = cache("E1", tag)
        add(5, f"e1_{name}", f"E1 ({name}): ten seeds, all-pairs gains, three rules, "
            "support resolution", sim("E1", src), [rev(src)], [src], eta="~1 h")
    src = cache("E6")
    add(6, "e6", "E6 word study: corrected held-out nulls, raw word discrepancies, "
        "independent Lemma S1 check", sim("E6", src), [rev(src)], [src], eta="~1 h")

    for key, tag, what in (("e4_base", BASE, "baseline prior: size + power"),
                           ("e4_prior", PRIOR, "matched prior: size only")):
        src = cache("E4", tag)
        add(7, key, f"E4 ({what}) with K-specific Test-3 strata; centres and "
            "projected-moment summaries kept", sim("E4", src), [rev(src)], [src],
            eta="~10 h" if key == "e4_base" else "~2 h")
    src = cache("E1")
    add(8, "selreps", "E1 selection on 10 x 10 evaluation corpora, three rules",
        sim("E1_selreps", src, "R_sel=10"),
        [f"Data/E1/e1_selreps_E1_full_{BASE}{sfx}.qs2"], [src], eta="~5 h")

    src = cache("E5")
    add(9, "e5_resolution", "support resolution of the E5 designs",
        ["Rscript", "Code/run_revision_diagnostics.R", str(w), "mode=e5",
         f"input={src}", f"out_suffix={sfx}"],
        [f"Data/E5/e5_resolution{sfx}.qs2"], [src], eta="~10 min")
    add(10, "unbinned", "protocol-matched UNBINNED comparator, phi floored at "
        "1e-14 / 1e-12 / 1e-10 and renormalised",
        ["Rscript", "Code/run_revision_diagnostics.R", str(w), "mode=unbinned",
         f"out_suffix={sfx}"],
        [f"Data/MDNA/unbinned_comparator{sfx}.qs2"],
        [PREP, cache("E1"), cache("E1", DGP2)], eta="~2.5 h")
    add(11, "restarts", "restart dispersion on the support harmonised over the "
        "production grid AND the restarts",
        ["Rscript", "Code/run_mdna_restarts.R", f"workers={w}", "common=2",
         f"out_suffix={sfx}"],
        [f"Data/MDNA/mdna_restarts_MDNA_2015_2016{sfx}.qs2"], [PREP], eta="~1.5 h",
        blas="default")

    for key, script, src in (("post_e4_base", "e4", cache("E4")),
                             ("post_e4_prior", "e4", cache("E4", PRIOR)),
                             ("post_e2_base", "e2_gap", cache("E2")),
                             ("post_e2_dgp2", "e2_gap", cache("E2", DGP2))):
        add(12, key, f"post-process {Path(rev(src)).name}",
            ["Rscript", f"Code/postprocess_{script}.R", f"file={rev(src)}"],
            [], [rev(src)], deps=[key.replace("post_", "")],
            code=[f"Code/postprocess_{script}.R"], eta="<1 min")
    production = [s["key"] for s in S if 2 <= s["n"] <= 12]
    add(13, "acceptance", "acceptance checks on every production output",
        ["Rscript", "Code/finalize_revision.R", f"suffix={sfx}"],
        [f"Results/revision{sfx}/acceptance.json"], [], deps=production,
        code=["Code/finalize_revision.R"], eta="~2 min")
    return S


def core_code():
    """Shared modules + configuration. A stage's identity adds the script(s) it
    runs (mirrors revision_code_files() in R), so a fix to one driver never
    invalidates a stage that ran another."""
    files = sorted([*Path("Code/R").glob("*.R"), Path("Code/config/configs.R")])
    return {str(p): sha256(p) for p in files}


def fits_state(manifest_csv):
    """sha256 of every cached fit now, checked against the pre-revision manifest."""
    now = {str(p): sha256(p) for p in sorted(Path("Data/FITS").glob("*.qs2"))}
    changed, missing = [], []
    if manifest_csv.exists():
        with open(manifest_csv, newline="") as f:
            for row in csv.DictReader(f):
                p = row["path"]
                if p not in now:
                    missing.append(p)
                elif now[p] != row["sha256"]:
                    changed.append(p)
    return now, changed, missing


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--suffix", default="_rev2")
    ap.add_argument("--select", default="all",
                    help="all | A | B | final | ids/keys 'a,b' | id range '2-6'")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    sfx = a.suffix
    if a.workers < 1 or not sfx.startswith("_") or not sfx.replace("_", "").isalnum():
        sys.exit("workers must be >= 1 and the suffix must look like _rev2")

    stages = build_stages(a.workers, sfx)
    by_key = {s["key"]: s for s in stages}
    sel = a.select
    if sel == "all":
        want = {s["key"] for s in stages}
    elif sel in BLOCKS:
        want = {s["key"] for s in stages if s["n"] in BLOCKS[sel]}
    elif "-" in sel and all(t.isdigit() for t in sel.split("-")):
        lo, hi = map(int, sel.split("-"))
        want = {s["key"] for s in stages if lo <= s["n"] <= hi}
    else:
        toks = sel.split(",")
        want = {s["key"] for s in stages if s["key"] in toks or str(s["n"]) in toks}
        unknown = [t for t in toks if t not in by_key and not any(str(s["n"]) == t for s in stages)]
        if unknown:
            sys.exit(f"unknown stage(s): {unknown}")
    if any(by_key[k]["n"] >= 2 for k in want):
        want |= {"units", "smoke"}           # gates are prerequisites, never optional

    dest = Path("Results") / ("revision" + sfx)
    dest.mkdir(parents=True, exist_ok=True)
    statefile = dest / "stages.json"
    state = json.loads(statefile.read_text()) if statefile.exists() else {}

    core = core_code()
    pkg = subprocess.check_output(
        ["Rscript", "--vanilla", "-e",
         'd <- packageDescription("OpTop"); cat(as.character(packageVersion("OpTop")), '
         'if (is.null(d$RemoteSha)) "" else d$RemoteSha)'], text=True).strip()

    def identity(s):
        cmd = [t for t in s["cmd"]]
        # the number of workers does not define a result
        cmd = ["<workers>" if t == str(a.workers) or t == str(min(a.workers, 2))
               or t == f"workers={a.workers}" else t for t in cmd]
        code = dict(core)
        scripts = [t for t in s["cmd"] if t.startswith("Code/") and t.endswith(".R")]
        code.update({p: sha256(p) for p in scripts + s["code"] if Path(p).is_file()})
        return dict(command=cmd, package=pkg, code=code,
                    inputs={p: sha256(p) for p in s["inputs"] if Path(p).is_file()})

    def valid(key):
        v = state.get(key, {})
        if v.get("status") != "complete":
            return False
        s = by_key[key]
        if any(not Path(p).is_file() for p in s["inputs"]):
            return False
        ident = identity(s)
        if any(v.get(f) != ident[f] for f in ("command", "package", "code", "inputs")):
            return False
        return all(Path(p).is_file() and sha256(p) == h
                   for p, h in v.get("outputs", {}).items()) and \
            all(p in v.get("outputs", {}) for p in s["outputs"])

    ordered = sorted(stages, key=lambda s: (ORDER.index(s["n"]), stages.index(s)))

    if a.dry_run:
        log(f"DRY RUN | suffix {sfx} | workers {a.workers} | OpTop {pkg}")
        for s in ordered:
            mark = "selected" if s["key"] in want else "-"
            log(f"  [{s['n']:>2}] {s['key']:<14} {('VALID' if valid(s['key']) else 'to run'):<7}"
                f" {mark:<9} {s['eta']:<8} {s['title']}")
        return 0

    lock = dest / "running.lock"
    if lock.exists():
        try:
            pid = int(lock.read_text().strip() or 0)
            os.kill(pid, 0)
            sys.exit(f"another batch is running (pid {pid}); lock: {lock}")
        except (ValueError, ProcessLookupError, PermissionError):
            log(f"removing stale lock {lock}")
            lock.unlink()
    lock.write_text(str(os.getpid()))
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))

    ran, verified, failed, blocked, skipped = [], [], [], [], []
    try:
        save_json(dest / "inputs.json", dict(
            suffix=sfx, package=pkg, core_code=core,
            outputs={s["key"]: s["outputs"] for s in stages},
            commands={s["key"]: s["cmd"] for s in stages}))
        manifest = Path("Results/csv/fit_manifest_rev1.csv")
        log("hashing the fit cache (before)")
        fits0, changed, missing = fits_state(manifest)
        if changed or missing:
            sys.exit(f"FIT CACHE DIFFERS FROM THE MANIFEST: {len(changed)} changed, "
                     f"{len(missing)} missing (first: {(changed + missing)[0]})")
        log(f"  {len(fits0)} fit files, all manifest entries unchanged")

        env0 = dict(os.environ)
        for s in ordered:
            k = s["key"]
            if k not in want:
                skipped.append(k)
                continue
            if valid(k):
                verified.append(k)
                log(f"VERIFIED {k} (complete under the current code, package and inputs)")
                continue
            gates = ["units", "smoke"] if s["n"] >= 2 else (["units"] if s["n"] == 1 else [])
            bad = [d for d in gates + s["deps"] if not valid(d)]
            if bad:
                blocked.append(k)
                state[k] = dict(status="blocked", blocked_by=bad,
                                at=datetime.now().isoformat(timespec="seconds"))
                save_json(statefile, state)
                log(f"BLOCKED  {k}: not complete -> {bad}")
                continue
            absent = [p for p in s["inputs"] if not Path(p).is_file()]
            if absent:
                failed.append(k)
                state[k] = dict(status="failed", reason=f"missing inputs {absent}")
                save_json(statefile, state)
                log(f"FAIL     {k}: missing inputs {absent}")
                continue
            logfile = dest / f"{k}.log"
            ident = identity(s)
            state[k] = dict(ident, status="running", log=str(logfile),
                            started=datetime.now().isoformat(timespec="seconds"))
            save_json(statefile, state)
            log(f"START    [{s['n']}] {k}: {s['title']} (expected {s['eta']}; log {logfile})")
            env = dict(env0)
            for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
                if s["blas"] == "default":
                    env.pop(v, None)
                else:
                    env[v] = s["blas"]
            if s["n"] >= 2:
                env["OPTOP_NO_FIT"] = "1"      # a cache miss is an error, never a fit
            else:
                env.pop("OPTOP_NO_FIT", None)  # gates fit toy models, as the unit suite does
            t0 = datetime.now()
            with logfile.open("a") as f:
                f.write(f"\n===== {t0.isoformat(timespec='seconds')} {' '.join(s['cmd'])}\n")
                f.flush()
                r = subprocess.run(s["cmd"], stdout=f, stderr=subprocess.STDOUT, env=env)
            outs_ok = all(Path(p).is_file() for p in s["outputs"])
            ok = r.returncode == 0 and outs_ok
            mins = (datetime.now() - t0).total_seconds() / 60
            state[k].update(
                status="complete" if ok else "failed", returncode=r.returncode,
                minutes=round(mins, 1),
                finished=datetime.now().isoformat(timespec="seconds"),
                outputs={p: sha256(p) for p in s["outputs"] if Path(p).is_file()})
            if r.returncode == 0 and not outs_ok:
                state[k]["reason"] = "exit 0 but a declared output is missing"
            save_json(statefile, state)
            (ran if ok else failed).append(k)
            log(f"{'DONE    ' if ok else 'FAIL    '} {k} in {mins:.1f} min"
                + ("" if ok else f" (exit {r.returncode}; see {logfile})"))
            if not ok and s["n"] <= 1:
                log("a gate failed: no production stage will run")
                break

        log("hashing the fit cache (after)")
        fits1, changed, missing = fits_state(manifest)
        drift = [p for p, h in fits0.items() if fits1.get(p) != h]
        new = sorted(set(fits1) - set(fits0))
        fits_ok = not (changed or missing or drift)
        production = [s["key"] for s in stages if s["n"] >= 2]
        summary = dict(
            suffix=sfx, selected=sorted(want), ran_now=ran, verified_complete=verified,
            failed=failed, blocked=blocked, skipped_not_selected=skipped,
            fits_unchanged=fits_ok, new_fit_files=new,
            production_complete=[k for k in production if valid(k)],
            production_outstanding=[k for k in production if not valid(k)],
            at=datetime.now().isoformat(timespec="seconds"))
        save_json(dest / "summary.json", summary)
        log(f"ran now          : {ran}")
        log(f"verified complete: {verified}")
        log(f"FAILED           : {failed}")
        log(f"BLOCKED          : {blocked}")
        log(f"skipped (not selected): {skipped}")
        log(f"fit cache unchanged: {fits_ok}"
            + (f" | new fit files (toy smoke fits only): {len(new)}" if new else ""))
        if not fits_ok:
            log("FIT IMMUTABILITY FAILURE -- results of this run must not be used")
        out = summary["production_outstanding"]
        if out:
            log(f"PRODUCTION STAGES STILL OUTSTANDING ({len(out)}): {out}")
        else:
            log("every production stage is complete and verified. This does NOT mean "
                "the manuscript is updated: exhibits, text and provenance come next.")
        return 1 if (failed or blocked or not fits_ok) else 0
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
