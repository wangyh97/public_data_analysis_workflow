"""Run independent samples with bounded processes and per-sample logs.
Usage: python scripts/23_dispatch_cellchat.py --workers 4
Existing completed result checkpoints are reused by script 21.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,json,subprocess,csv,re
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);args=p.parse_args()
assert 1<=args.workers<=6
cfg=json.loads((ROOT/'analysis_config.json').read_text())
jobs=list(csv.DictReader((ROOT/'data/CellChat/jobs.csv').open()))
def run(j):
 job=j['job'];log=ROOT/'audit'/('21_dispatch_'+job+'.log')
 with log.open('w') as f:
  r=subprocess.run([cfg['paths']['Rscript'],str(ROOT/'scripts/21_run_cellchat.R'),'^'+re.escape(job)+'$'],stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError(f'{job} failed; see {log}')
 return job
with ThreadPoolExecutor(max_workers=args.workers) as ex:
 for job in ex.map(run,jobs):print('COMPLETE',job,flush=True)
