def timestamp_fn():
    return datetime.now().strftime("%Y_%m_%d_%H_%M_%S")


def ensure_dir_fn(path: Path, verbose: bool = False):
    if path.exists():
        if verbose:
            print(f"[INFO] Output directory already exists: {path}")
    else:
        if verbose:
            print(f"[INFO] Creating output directory: {path}")
        path.mkdir(parents=True, exist_ok=True)


def open_logfile_fn(outdir: Path):
    logfile = outdir / "amplicon_seq.log"
    return logfile.open("w")


def write_log(fh, msg: str, verbose=True):
    print(msg, file=fh)
    if verbose:
        print(msg, file=sys.stderr)