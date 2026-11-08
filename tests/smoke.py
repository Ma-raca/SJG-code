#!/usr/bin/env python3
"""Build/search integration and independent distance/ground-truth checks."""
import importlib.util, math, os, subprocess, sys, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('generate', ROOT/'examples/generate.py')
generator = importlib.util.module_from_spec(spec); spec.loader.exec_module(generator)
binary = Path(sys.argv[1] if len(sys.argv)>1 else ROOT/'sjg').resolve()
def edit(a,b):
    prev = list(range(len(b)+1))
    for i, x in enumerate(a,1):
        cur = [i]
        for j,y in enumerate(b,1): cur.append(min(cur[-1]+1, prev[j]+1, prev[j-1]+(x!=y)))
        prev=cur
    return prev[-1]
def distance(a,b,w):
    v=math.sqrt(sum((int(x)-int(y))**2 for x,y in zip(a[0].split(','),b[0].split(','))))/50
    x,y=set(a[3].split('|')),set(b[3].split('|'))
    return sum(p*q for p,q in zip(w,[v,abs(int(a[1])-int(b[1]))/50,edit(a[2],b[2])/6,1-len(x&y)/len(x|y)]))
def answers(p):
    rows=[]
    for line in p.read_text().splitlines():
        s=line.split()
        if s[0]=='k': rows.append([int(x) for x in s[3:]])
        elif s[0]=='D':
            assert len(s)-1==len(rows[-1])
    return rows
with tempfile.TemporaryDirectory(prefix='sjg-test-') as tmp:
    d=Path(tmp); generator.generate(d)
    env=dict(os.environ,OMP_NUM_THREADS='1',HJG_THREADS='1',HJG_SKYLINE='1',SKY_FLOOR='16',SKY_MODE='skyline',SJG_INDEX_K='15',ND_K='16',ND_EMIT='15',ND_ITERS='12',OUT_DIST='1')
    for key in ['EF_LIST','EXTERNAL_QUERY_DATA','NONORM','NOMYERS','SKY_CAP','SKY_ALPHA','SKY_ALPHA_M','NO_NATIVE_OUT']:env.pop(key,None)
    def run(args): subprocess.run([str(binary)]+[str(x) for x in args],env=env,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    run(['build',d/'database.txt',d/'graph.ind',d/'graph.sjg'])
    db=[s.split() for s in (d/'database.txt').read_text().splitlines()[3:]]
    qs=[s.split() for s in (d/'query_objects.txt').read_text().splitlines()[3:]]
    env['EXTERNAL_QUERY_DATA']=str(d/'query_objects.txt')
    for w in [(0.25,0.25,0.25,0.25),(1,0,0,0),(0.1,0.2,0.3,0.4)]:
        (d/'queries.txt').write_text(' '.join(map(str,w))+'\n1\n10\n0\n16\n'+' '.join(map(str,range(16)))+'\n')
        env['OUT']=str(d/'truth.txt');run(['brute',d/'database.txt',d/'queries.txt'])
        exact=answers(d/'truth.txt');assert len(exact)==16
        for q,ids in zip(qs,exact):
            assert len(set(ids))==10 and all(0<=i<len(db) for i in ids)
            threshold=sorted(distance(q,x,w) for x in db)[9]
            assert all(distance(q,db[i],w)<=threshold+1e-6 for i in ids)
        env['OUT']=str(d/'answers.txt');env['EF']='256';run(['search',d/'database.txt',d/'graph.sjg',d/'queries.txt'])
        found=answers(d/'answers.txt');assert len(found)==16
        for ids in found: assert len(set(ids))==10 and all(0<=i<len(db) for i in ids)
        # Accept tied top-k memberships using independently computed distances.
        hits=0
        for q,ids in zip(qs,found):
            threshold=sorted(distance(q,x,w) for x in db)[9]
            hits+=sum(distance(q,db[i],w)<=threshold+1e-6 for i in ids)
        assert hits/160>=.95, ('unexpectedly low smoke recall',w,hits/160)
    print('PASS: build, external queries, exact scan, three weight profiles, independent distance oracle.')
