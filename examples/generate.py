#!/usr/bin/env python3
"""Generate a reproducible, synthetic four-component dataset (no downloads)."""
import argparse, random
from pathlib import Path

def generate(root, n=256, nq=16):
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261108)
    def row():
        vector = ','.join(str(rng.randrange(30)) for _ in range(3))
        scalar = str(rng.randrange(50))
        name = ''.join(rng.choice('abcdef') for _ in range(6))
        tags = '|'.join(sorted(rng.sample(['a','b','c','d','e','f'], 3)))
        return f'{vector} {scalar} {name} {tags}\n'
    meta = '2 1 6 7\n50 50 6 1\n'
    (root/'database.txt').write_text(f'{n} 4\n'+meta+''.join(row() for _ in range(n)))
    (root/'query_objects.txt').write_text(f'{nq} 4\n'+meta+''.join(row() for _ in range(nq)))
    (root/'queries.txt').write_text('0.25 0.25 0.25 0.25\n1\n10\n0\n'+str(nq)+'\n'+' '.join(map(str, range(nq)))+'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', default='output/example'); args = parser.parse_args()
    generate(args.output)
    print('Generated database, external query objects, and query protocol in', args.output)
