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

Everything below funnels into [`run_grid.sh`](run_grid.sh), which builds the
engine and writes `data/N<N>/p<p8>/data.dat` — the exact layout
`tools/analyze.py --all` and `tools/heatmap.py` expect. Do the analysis/plots
locally after downloading `data/`.

## Option A — AWS EC2 spot (recommended one-off, ~$3–5)

A 96-vCPU compute instance for ~1–2 h grinds the whole grid. `c7a.24xlarge`
(AMD Zen4, x86) or `c8g.24xlarge` (Graviton4, ARM — usually cheaper; the C code
is portable, just builds with gcc). Prerequisite: AWS CLI configured
(`aws configure`) and an SSH key pair.

```sh
# 1. latest Ubuntu 24.04 AMI for your region (x86; use .../arm64/... for c8g)
AMI=$(aws ssm get-parameters --names \
  /aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id \
  --query 'Parameters[0].Value' --output text)

# 2. request a spot instance
aws ec2 run-instances \
  --image-id "$AMI" --instance-type c7a.24xlarge \
  --instance-market-options 'MarketType=spot' \
  --key-name YOUR_KEY --security-groups YOUR_SG \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=brane}]'

# 3. rsync the repo up, run, pull results back  (IP from the console/CLI)
IP=<public-ip>
rsync -az --exclude data --exclude .git ./ ubuntu@$IP:brane/
ssh ubuntu@$IP 'cd brane && NT=96 SWEEPS=4000 EPS=0.005 bash cloud/run_grid.sh'
rsync -az ubuntu@$IP:brane/data/ ./data/

# 4. TERMINATE (spot instances still cost while running!)
aws ec2 terminate-instances --instance-ids <id>
```

Fully-unattended variant: pass `run_grid.sh` as `--user-data` with
`S3_BUCKET=your-bucket SHUTDOWN=1`, and it uploads `brane_data_*.tgz` to S3 and
powers off when done. (Needs an instance profile with `s3:PutObject`.)

## Option B — GCP spot

```sh
gcloud compute instances create brane \
  --machine-type c3-highcpu-176 --provisioning-model=SPOT \
  --image-family ubuntu-2404-lts-amd64 --image-project ubuntu-os-cloud
gcloud compute scp --recurse ./ brane:~/brane
gcloud compute ssh brane --command 'cd brane && bash cloud/run_grid.sh'
gcloud compute scp --recurse brane:~/brane/data ./data
gcloud compute instances delete brane        # don't forget
```

## Option C — university HPC (SLURM) — free, do this if you have access

[`slurm_grid.sbatch`](slurm_grid.sbatch) is an array job: one task per `(N, p8)`
cell, each using its allocated cores as replicas.

## Option D — Apple Simcloud (ACS) — internal batch compute

See **[`SIMCLOUD.md`](SIMCLOUD.md)** for the full guide. One job per `(N, p8)`
cell via `simcloud batch post`; use the `mr2-as` M2 Ultra cluster (~3.6× faster
per core than x86, near-linear scaling). Driver scripts: `simcloud_submit.sh`,
`simcloud_fetch.sh`, `simcloud_bench_run.sh`.

```sh
# edit NS / P8S and set --array=0-(len(NS)*len(P8S)-1), then:
sbatch cloud/slurm_grid.sbatch
# when finished, back on your laptop:
uv run tools/analyze.py --all
uv run tools/heatmap.py --replot-all --png plots/heatmap_combined.png
uv run tools/analyze.py --combined
```

## Tuning `run_grid.sh`

All via environment variables (see the top of the script): `NS`, `P8S`, `NT`
(default = all cores), `THERM`, `SWEEPS`, `EPS` (convergence target on the
`Delta2` relative error — lower = tighter = longer), `MINSW`, `S3_BUCKET` /
`GCS_BUCKET` for auto-upload, `SHUTDOWN=1` to poweroff on completion.

For the **statistics** goal, keep `NT` = all cores and lower `EPS` (e.g. 0.005).
For the **large-N** goal, add bigger sizes to `NS` and prefer a high-clock box.
