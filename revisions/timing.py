"""Time complete fits; keep worker latency separate from batch throughput."""
from .common import *
import argparse
import time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor


_WARMUP = None

def initialize(warmup_task=None):
    # Import all comparators before measurement; warm-up is separately recorded.
    from peak_extraction import extract_peaks
    global _WARMUP
    if warmup_task is not None:
        _WARMUP = timed(warmup_task)


def timed(task):
    method, row_id, x, y = task
    start = time.perf_counter_ns()
    try:
        peak = fit(method, x, y)
        state, error = ('ok', '') if np.isfinite(peak) else ('fit_error', 'Nonfinite peak')
    except Exception as exc:
        peak, state, error = np.nan, status(exc), str(exc)
    return dict(method=method, source_row=row_id, worker_pid=os.getpid(), warmup_ms=None if _WARMUP is None else _WARMUP["elapsed_ms"], elapsed_ms=(time.perf_counter_ns()-start)/1e6,
                peak=peak, status=state, error=error)


def summarize(directory):
    directory = Path(directory)
    raw = pd.read_csv(directory/'timings.csv')
    batches = pd.read_csv(directory/'batches.csv')
    summary = raw.groupby(['method','workers']).elapsed_ms.agg(
        count='size', median_ms='median', p95_ms=lambda a:a.quantile(.95),
        p99_ms=lambda a:a.quantile(.99), maximum_ms='max').reset_index()
    summary = summary.merge(batches.groupby(['method','workers']).effective_ms_per_trace.agg(
        effective_ms_per_trace='median').reset_index())
    summary.to_csv(directory/'timing_table.csv', index=False)
    print(summary.to_string(index=False))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--workers', type=int, nargs='+', default=[1,18])
    p.add_argument('--methods', nargs='+', choices=METHODS, default=METHODS)
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--hardware-label', default=platform.node())
    p.add_argument('--limit', type=int, help='Smoke test only; omit for manuscript timings')
    p.add_argument('--summarize-only', action='store_true')
    a=p.parse_args()
    if a.summarize_only: return summarize(a.out)
    if not a.input: p.error('--input is required')
    if min(a.workers)<1 or a.repeats<1: p.error('workers/repeats must be positive')
    a.out.mkdir(parents=True, exist_ok=True)
    df=structured(a.input)
    total=len(df)
    if a.limit: df=df.iloc[:a.limit]
    meta=provenance([a.input]);meta.update(hardware_label=a.hardware_label, total_input_traces=total,
        measured_traces=len(df), smoke_test=bool(a.limit), repeats=a.repeats,
        protocol='CSV parsing and imports excluded. Full fit including selection timed. Batch wall time includes startup, warmup, IPC and collection. Warmup latencies recorded separately; no worker-latency division by workers.')
    save_json(a.out/'metadata.json',meta)
    initialize(); records=[]; batches=[]; warmup=[]
    for workers in a.workers:
        for method in a.methods:
            tasks=[(method,r.source_row,r.voltage,r.current) for r in df.itertuples()]
            for repeat in range(a.repeats):
                start=time.perf_counter()
                if workers==1:
                    warmup.append(timed(tasks[0]))
                    rows=list(map(timed,tasks))
                else:
                    with ProcessPoolExecutor(workers, mp_context=mp.get_context('spawn'), initializer=initialize, initargs=(tasks[0],)) as pool:
                        rows=list(pool.map(timed,tasks,chunksize=1))
                wall=time.perf_counter()-start
                for r in rows: r.update(workers=workers,repeat=repeat)
                records.extend(rows)
                batches.append(dict(method=method,workers=workers,repeat=repeat,n=len(rows),wall_seconds=wall,effective_ms_per_trace=wall*1000/len(rows)))
                pd.DataFrame(records).to_csv(a.out/'timings.csv',index=False)
                pd.DataFrame(batches).to_csv(a.out/'batches.csv',index=False)
                print(method,workers,repeat,wall,flush=True)
    pd.DataFrame(warmup).to_csv(a.out/'warmup.csv',index=False)
    pd.DataFrame(records)[['method','workers','repeat','worker_pid','warmup_ms']].drop_duplicates().to_csv(a.out/'worker_warmup.csv',index=False)
    summarize(a.out)

if __name__=='__main__': main()
