import datetime

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

configfile: "config/config.yml"

OUTPUT_DIR = config["output_dir"]
FQLIST = config["fqlist"]
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
