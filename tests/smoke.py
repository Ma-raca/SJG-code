#!/usr/bin/env python3
"""Build/search integration and independent distance/ground-truth checks."""
import importlib.util, math, os, struct, subprocess, sys, tempfile
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
    for key in ['EF_LIST','EXTERNAL_QUERY_DATA','NONORM','NOMYERS','SKY_CAP','SKY_ALPHA','SKY_ALPHA_M','NO_NATIVE_OUT','SJG_EDGE_BOUNDS']:env.pop(key,None)
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
        for enabled in [False,True]:
            if enabled: env['SJG_EDGE_BOUNDS']='1'
            else: env.pop('SJG_EDGE_BOUNDS',None)
            env['OUT']=str(d/'answers.txt');env['EF']='256';run(['search',d/'database.txt',d/'graph.sjg',d/'queries.txt'])
            found=answers(d/'answers.txt');assert len(found)==16
            for ids in found: assert len(set(ids))==10 and all(0<=i<len(db) for i in ids)
            # Accept tied top-k memberships using independently computed distances.
            hits=0
            for q,ids in zip(qs,found):
                threshold=sorted(distance(q,x,w) for x in db)[9]
                hits+=sum(distance(q,db[i],w)<=threshold+1e-6 for i in ids)
            assert hits/160>=.95, ('unexpectedly low smoke recall',w,enabled,hits/160)

    # A one-dimensional star gives a known nearest neighbor without distance ties.
    # This also checks that the switch changes actual distance evaluations.
    (d/'line.txt').write_text('8 1\n1\n1\n'+'\n'.join(map(str,range(8)))+'\n')
    (d/'line_query.txt').write_text('1 1\n1\n1\n0.25\n')
    (d/'line_protocol.txt').write_text('1\n1\n1\n0\n1\n0\n')
    (d/'star.sjg').write_bytes(b'SJG1'+struct.pack('<4i',8,1,1,0)
        +struct.pack('<9q',0,*([7]*8))+struct.pack('<7i',*range(1,8))
        +struct.pack('<7f',*range(1,8)))
    env['EXTERNAL_QUERY_DATA']=str(d/'line_query.txt');env['EF']='2'
    counts={}
    for setting in [None,'0','1']:
        if setting is None: env.pop('SJG_EDGE_BOUNDS',None)
        else: env['SJG_EDGE_BOUNDS']=setting
        env['OUT']=str(d/'line_answers.txt')
        run(['search',d/'line.txt',d/'star.sjg',d/'line_protocol.txt'])
        assert answers(d/'line_answers.txt')==[[0]]
        counts[setting]=next(float(l.split()[1]) for l in (d/'line_answers.txt').read_text().splitlines() if l.startswith('dc'))
    assert counts[None]==counts['0'] and counts['1']<counts['0'],counts
    print('PASS: mixed-metric oracle; edge filtering defaults off; opt-in preserves the known neighbor and saves distance evaluations.')
