# Running brane on a cluster / cloud

The engine is **replica-parallel and CPU-bound**: one Markov chain per core,
averaged. So the statistical error on every `G(q)` point falls as
`~ 1/sqrt(nt * sweeps)`. Our laptop bottleneck is *samples*, not code — a 96-core
box gives `sqrt(96/12) ~ 2.8x` smaller error bars at the same wall-time, and the
scatter that keeps our fits off a clean line (see the main README / thesis
comparison) largely goes away.

Big-N *reach* is a different limit: a single chain's sweep costs `~N^4` and is
sequential, so it is bounded by single-core **clock**, not core count. Pick a
high-clock instance if you want to push N; pick a high-core-count instance if you
want statistics. Spot/preemptible pricing is ideal here — replicas are
independent, so an interrupted run is just discarded.

## Apple Simcloud (ACS) — the workflow we use

See **[`SIMCLOUD.md`](SIMCLOUD.md)** for the full guide (cluster comparison, VPC
net ids, sizing). One job per `(N, p8)` cell via `simcloud batch post`, on the
`mr2-as` M2 Ultra cluster (~3.6× faster per core than x86, near-linear scaling).

The whole loop is three steps — **launch → monitor → fetch**:

```sh
# 1. launch the production grid (edit knobs inline; see overnight.sh header)
bash cloud/overnight.sh
#    -> prints the exact monitor/fetch commands with the batch id filled in

# 2. watch it (live progress bar; uses cloud/.last_batch)
CLUSTER=mr2-as bash cloud/simcloud_monitor.sh
#    (or the built-in: simcloud -c mr2-as job wait --batch <id> --summary)

# 3. pull results into data/ once jobs finish
CLUSTER=mr2-as bash cloud/simcloud_fetch.sh

# 4. analyze locally
uv run tools/analyze.py --all
uv run tools/heatmap.py --replot-all --png plots/heatmap_combined.png
```

### Scripts

| script | role |
|---|---|
| `overnight.sh` | **entrypoint** — sets production defaults, then execs `simcloud_submit.sh` |
| `simcloud_submit.sh` | builds the source + toolchain bundles and posts the batch |
| `simcloud_task.sh` | runs **inside** each job: maps the batch index → one `(N,p8)` cell, builds, runs the engine |
| `simcloud_monitor.sh` | live progress-bar monitor of a batch |
| `simcloud_fetch.sh` | waits for the batch, downloads the output bundles, merges into `data/` |

> Recipes for non-Apple clusters (AWS/GCP spot, SLURM array jobs via
> `run_grid.sh` / `slurm_grid.sbatch`) were removed to keep this focused on the
> Simcloud path; they remain in the git history if ever needed again.
