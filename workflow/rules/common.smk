import datetime

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

configfile: "config/config.yml"

OUTPUT_DIR = config["output_dir"]
FQLIST = config["fqlist"]
PLATFORM=config["platform"]
GATK_BUNDLE_DIR = config["gatk_bundle_dir"]
REF = config["ref"]
REF_INDEX_PREFIX = config["ref_index_prefix"]

common_params = {
    "ref": REF,
    "gatk_bundle_dir": GATK_BUNDLE_DIR
}

wildcard_constraints:
    sample_id="CL1\d+L\d+",
    cid="CL1\d\d\d\d",
    lid="L\d\d",
    snn="\d\d"

def get_log_path(sample_id):
    return os.path.join(config["log_dir"], timestamp, f"{sample_id}.log")

# Get sample infor from falist
def get_samples(fqlist_path):
    samples = []
    with open(fqlist_path, "r") as f:
        for line in f:
            fq, cid, lid, samid = line.strip().split()
            samples.append((fq, cid, lid, samid))
    return samples

SAMPLES = get_samples(FQLIST)

def generate_regional_segments(chroms):
    """
    Generate chromosome segments as (chr_id, start, end) tuples.

    This function divides each chromosome into fixed-size regions (delta) based on the reference genome's .fai file
    and yields each segment as a tuple of (chr_id, start, end).

    Parameters:
        chroms: A list of chromosome IDs (e.g., ["chr1", "chr2"]) or a single chromosome ID (e.g., "chr1").

    Returns:
        A generator that yields tuples in the format (chr_id, start, end).
    """
    delta = config["chunk_size"]
    in_fai = config["ref_fai"]

    # Ensure chroms is a list (even if a single chromosome ID is provided)
    if isinstance(chroms, str):
        chroms = [chroms]

    with open(in_fai) as fh:
        for r in fh:
            col = r.strip().split()
            chr_id = col[0]
            chr_length = int(col[1])

            # Skip chromosomes not in the target list
            if chroms is not None and len(chroms):
                if chr_id not in chroms:
                    continue

            # Generate regions for the chromosome
            for i in range(0, chr_length, delta):
                start = i + 1
                end = i + delta if i + delta <= chr_length else chr_length
                yield (chr_id, start, end)
