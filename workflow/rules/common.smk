import json, datetime, pathlib

configfile: "config/config.yml"

stamp = pathlib.Path(".run_ts.json")
if stamp.exists():
    RUN_TS = json.loads(stamp.read_text())["run_ts"]
else:
    RUN_TS = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    stamp.write_text(json.dumps({"run_ts": RUN_TS}))

onsuccess:
    stamp.unlink(missing_ok=True)

onerror:
    stamp.unlink(missing_ok=True)

config["run_ts"] = RUN_TS

OUTPUT_DIR = config["output_dir"]
FQLIST = config["sample_list"]
PLATFORM=config["platform"]
GATK_BUNDLE_DIR = config["gatk_bundle_dir"]
REF = config["ref"]
BENCH_DIR = config["benchmark_dir"]

common_params = {
    "ref": REF,
    "gatk_bundle_dir": GATK_BUNDLE_DIR
}

wildcard_constraints:
    start="\d+",
    end="\d+",
    chr_id="chr.+",
    chrom="chr.+"

def get_log_path(sample_id):
    return os.path.join(config["log_dir"], config["run_ts"], f"{sample_id}.log")

def get_samples(fqlist_path):
    samples = []
    with open(fqlist_path, "r") as f:
        for line in f:
            stripped_line = line.strip()
            if stripped_line:
                samples.append(stripped_line)
    return samples

SAMPLES = get_samples(FQLIST)

def generate_chromosome_segments(exclude=True, chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]):
    """
    Generate chromosome segments as (chr_id, start, end) tuples, excluding specified regions.

    This function divides each chromosome into fixed-size regions (delta) based on the reference genome's .fai file.
    If exclude_regions are provided, the function splits the chromosome segments to avoid overlapping with excluded regions,
    ensuring full coverage of the remaining regions.

    exclude_regions: A list of regions to exclude, in the format ["chr:start-end", ...].
                    Example: ["chr21:1-5000000", "chr22:1-5000000"].

    Returns:
        A generator that yields tuples in the format (chr_id, start, end).
    """
    delta = config["chunk_size"]
    in_fai = config["ref_fai"]
    exclude_regions = config["exclude_regions"] if exclude else None

    # Parse exclude_regions into a dictionary for quick lookup
    exclude_dict = {}
    if exclude_regions:
        for region in exclude_regions:
            chr_id, range_str = region.split(":")
            start, end = map(int, range_str.split("-"))
            if chr_id not in exclude_dict:
                exclude_dict[chr_id] = []
            exclude_dict[chr_id].append((start, end))

    with open(in_fai) as fh:
        for r in fh:
            col = r.strip().split()
            chr_id = col[0]
            chr_length = int(col[1])

            # Skip chromosomes not in the target list
            if chroms is not None and len(chroms):
                if chr_id not in chroms:
                    continue

            # Get excluded regions for this chromosome (if any)
            excluded_segments = exclude_dict.get(chr_id, [])

            # Sort excluded regions by start position
            excluded_segments.sort()

            # Initialize the start of the chromosome
            current_start = 1

            # Iterate over the chromosome, splitting around excluded regions
            while current_start <= chr_length:
                # Find the next excluded region that overlaps with the current segment
                next_excluded = None
                for excl_start, excl_end in excluded_segments:
                    if excl_start <= current_start <= excl_end or current_start <= excl_start <= current_start + delta:
                        next_excluded = (excl_start, excl_end)
                        break

                if next_excluded:
                    excl_start, excl_end = next_excluded

                    # Yield the segment before the excluded region
                    if current_start < excl_start:
                        yield (chr_id, current_start, excl_start - 1)

                    # Move current_start to the end of the excluded region
                    current_start = excl_end + 1
                else:
                    # No more excluded regions, yield the remaining segment
                    end = min(current_start + delta - 1, chr_length)
                    yield (chr_id, current_start, end)
                    current_start = end + 1

