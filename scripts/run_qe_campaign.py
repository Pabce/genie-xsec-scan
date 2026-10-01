#!/usr/bin/env python3
"""Checkpointed Ar/Ca QE scans; reuse completed chunks only with matching runtime."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import pathlib
import subprocess
import time
import fcntl

TARGETS = {'Ar40': 1000180400, 'Ca40': 1000200400, 'Ca48': 1000200480}


def sha(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def save_atomic(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scanner', type=pathlib.Path, required=True)
    parser.add_argument('--settings', type=pathlib.Path, required=True,
                        help='JSON list: E_MeV, theta_deg, omega_max_MeV')
    parser.add_argument('--output', type=pathlib.Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--chunk-size', type=int, default=25)
    parser.add_argument('--tolerances', nargs='+', type=float, default=[.001, .00025])
    args = parser.parse_args()
    if args.jobs < 1 or args.chunk_size < 1:
        parser.error('jobs and chunk-size must be positive')
    if any(not 0 < t < 1 for t in args.tolerances):
        parser.error('tolerances must be finite and between zero and one')
    args.scanner = args.scanner.resolve(strict=True)
    args.settings = args.settings.resolve(strict=True)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / '.campaign.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    genie = pathlib.Path(os.environ['GENIE']).resolve(strict=True)
    roots = {genie / 'config', args.scanner.parent.parent / 'config'}
    roots.update(pathlib.Path(x).resolve() for x in os.environ.get('GXMLPATH', '').split(':') if x)
    files = {args.scanner, args.settings}
    for root in roots:
        files.update(p.resolve() for p in root.rglob('*.xml'))
    files.update(p.resolve() for p in (genie / 'lib').glob('*.dylib'))
    files.update(p.resolve() for p in (genie / 'lib').glob('*.so*'))
    identity = {'files': {str(p): sha(p) for p in sorted(files)},
                'GENIE': str(genie), 'GXMLPATH': os.environ.get('GXMLPATH', '')}
    identity_sha = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    manifest = args.output / 'runtime.json'
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise RuntimeError('Runtime/settings changed; use a new versioned output directory')
    else:
        save_atomic(manifest, identity)
    settings = json.loads(args.settings.read_text())
    tasks = []
    for s in settings:
        e, angle, end = s['E_MeV'], s['theta_deg'], s['omega_max_MeV']
        if not (math.isfinite(e) and e > end >= 1 and isinstance(end, int)
                and math.isfinite(angle) and 0 < angle < 180):
            raise ValueError(f'Invalid setting: {s}')
        for nucleus in TARGETS:
            for tolerance in args.tolerances:
                for start in range(1, end + 1, args.chunk_size):
                    tasks.append((nucleus, e, angle, start,
                                  min(start + args.chunk_size - 1, end), tolerance))
    if len(set(tasks)) != len(tasks):
        raise ValueError('Duplicate settings/tolerances would create duplicate work')

    def run(task):
        nucleus, e, angle, start, end, tolerance = task
        tag = f'{nucleus}_E{e:g}_th{angle:g}_w{start}-{end}_tol{tolerance:g}'
        output = args.output / (tag + '.csv')
        receipt = args.output / (tag + '.receipt.json')
        partial = args.output / (tag + '.partial.csv')
        log = args.output / (tag + '.log')
        command = [str(args.scanner), '--mode', 'tune', '--tune', 'G18_10a_02_11a',
                   '--event-generator-list', 'EMQE', '--probe', '11',
                   '--target', str(TARGETS[nucleus]), '--observable', 'd2',
                   '--diff', 'Eprime,costheta_l', '--fixed', f'E={e / 1000:.17g}',
                   '--fixed', f'costheta_l={math.cos(math.radians(angle)):.17g}',
                   '--scan', f'Eprime:{(e-start)/1000:.17g}:{(e-end)/1000:.17g}:{end-start+1}',
                   '--xsec-unit', 'nb', '--components', '--fold', 'auto', '--strict',
                   '--qel-fold-density', 'adaptive-theta', '--qel-rel-tol', str(tolerance),
                   '--output', str(partial)]
        if receipt.exists():
            previous = json.loads(receipt.read_text())
            if (previous['runtime_sha256'] != identity_sha or previous['command'] != command
                    or previous['csv_sha256'] != sha(output)):
                raise RuntimeError(f'Changed completed artifact: {tag}')
            return tag, 'reused'
        if output.exists() or partial.exists() or log.exists():
            raise RuntimeError(f'Unreceipted artifact preserved; inspect before retrying: {tag}')
        begin = time.monotonic()
        with log.open('x') as stream:
            process = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        if process.returncode:
            raise RuntimeError(f'Failed {tag}; partial output and log preserved')
        with partial.open() as stream:
            rows = list(csv.DictReader(line for line in stream if not line.startswith('#')))
        totals = [row for row in rows if row['component'] == 'total']
        if (len(totals) != end-start+1 or
                any(row['status'] != 'ok' or not math.isfinite(float(row['value']))
                    or float(row['value']) < 0 for row in rows)):
            raise RuntimeError(f'Invalid or incomplete output: {tag}')
        if sorted(round(float(row['omega']) * 1000) for row in totals) != list(range(start, end+1)):
            raise RuntimeError(f'Incorrect output grid: {tag}')
        partial.replace(output)
        save_atomic(receipt, {'runtime_sha256': identity_sha, 'command': command,
                             'csv_sha256': sha(output), 'points': len(totals),
                             'elapsed_seconds': time.monotonic()-begin})
        return tag, 'completed'

    # Bound the submission queue. On a failure, finish/checkpoint already-running
    # chunks, but launch no further work. No failed scientific result is retried.
    todo = iter(tasks)
    begin, count, failures = time.monotonic(), 0, []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {}
        def submit():
            task = next(todo, None)
            if task is not None:
                pending[pool.submit(run, task)] = task
        for _ in range(args.jobs):
            submit()
        while pending:
            ready, _ = concurrent.futures.wait(pending, timeout=20,
                                               return_when=concurrent.futures.FIRST_COMPLETED)
            if not ready:
                print(f'Progress: {count}/{len(tasks)} chunks; {len(pending)} active; '
                      f'{time.monotonic()-begin:.0f}s elapsed', flush=True)
            for future in ready:
                del pending[future]
                try:
                    tag, state = future.result()
                    count += 1
                    print(f'{state}: {tag} ({count}/{len(tasks)})', flush=True)
                except Exception as error:
                    failures.append(str(error))
                    print(f'ERROR: {error}', flush=True)
            if not failures:
                for _ in ready:
                    submit()
    if failures:
        raise RuntimeError('\n'.join(failures))
    print(f'Complete: {count} verified chunks', flush=True)


if __name__ == '__main__':
    main()
