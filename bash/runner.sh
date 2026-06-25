#!/bin/bash
# == Information ==============================================================
#  Run the pipeline for processing the data
#  Change the code and loop as the scripts required


# == Configuration ============================================================
set -euo pipefail


# == Runner ===================================================================
# Name of the samples
SAMPLES=(siC_1 siC_2 siC_3 siBRG1_1 siBRG1_2 siUAP56_1 siUAP56_2)

# Running loop
for SAMPLE in "${SAMPLES[@]}"; do
    for STRAND in fwd rev; do
        bash 08_peaks.sh "${SAMPLE}" "${STRAND}"
    done
done

