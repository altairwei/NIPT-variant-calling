SCRIPT_DIR="$(cd "$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")" && pwd)" || die "Couldn't determine the script's running directory, which probably matters, bailing out" 2

export PATH="${SCRIPT_DIR}/bin:${SCRIPT_DIR}/scripts:$PATH"

source activate nipt
