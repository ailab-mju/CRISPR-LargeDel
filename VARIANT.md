# cli/ variant: detection

A copy of `cli/` (v20261023) with these switches added, **on by default**: --ld_confirm (#25, default 1), --split_wt_margin (#26, default 2).

Turn any of them off on the command line to get `cli/`'s behaviour for that part
(`--ld_confirm 0`, `--split_wt_margin -1`).
With all of them off the output is identical to `cli/`.

The code is the one tested in `comparison_results/cores/cli_detect`; results and side effects are in
`comparison_results/detect_tests/README.md` and the ExCas Core Changes summary.
Made by `comparison_results/scripts/make_cli_variant.py detection`.

Build: `make HTS=0` (native, writes SAM) or `make` (with htslib, BAM).
Web core: `wasm/build_wasm.sh` (emsdk 3.1.18) writes `wasm/out/kmer_align.js` + `kmer_align.wasm`,
the files the website's `docs/` folder loads.
