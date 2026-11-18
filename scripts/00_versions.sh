#!/usr/bin/env bash
# Print the versions of every tool used, for the Methods / README.
source "$(dirname "$0")/../config/config.sh"
echo "GATK:     $($GATK --version 2>/dev/null | grep -m1 -i 'gatk')"
echo "bcftools: $(bcftools --version | head -1)"
echo "samtools: $(samtools --version | head -1)"
echo "tabix:    $(tabix --version 2>&1 | head -1)"
echo "WhatsHap: $(whatshap --version 2>/dev/null)"
echo "Python:   $(python3 --version)"
python3 -c "import pyBigWig, matplotlib; print('pyBigWig', pyBigWig.__version__, '| matplotlib', matplotlib.__version__)"
echo "VEP image: $VEP_IMAGE"
