#!/usr/bin/bash

job_id=$(echo $1 | cut -f1 -d ";")

job_state_code=$(scontrol show job $job_id | \
                 grep -oh "JobState=[[:alpha:]]*" | \
                 cut -d "=" -f2)

# job_state_code=$(squeue -j $job_id -o "%t" | \
#                  tail -n 1)

declare -A code_to_state=( ["BOOT_FAIL"]="failed" \
                           ["CANCELLED"]="failed" \
                           ["COMPLETED"]="success" \
                           ["CONFIGURING"]="running" \
                           ["COMPLETING"]="running" \
                           ["DEADLINE"]="failed" \
                           ["FAILED"]="failed" \
                           ["NODE_FAIL"]="failed" \
                           ["OUT_OF_MEMORY"]="failed" \
                           ["PENDING"]="running" \
                           ["PREEMPTED"]="failed" \
                           ["RUNNING"]="running" \
                           ["RESV_DEL_HOLD"]="running" \
                           ["REQUEUE_FED"]="running" \
                           ["REQUEUE_HOLD"]="running" \
                           ["REQUEUED"]="running" \
                           ["RESIZING"]="running" \
                           ["REVOKED"]="failed" \
                           ["SIGNALING"]="running" \
                           ["SPECIAL_EXIT"]="failed" \
                           ["STAGE_OUT"]="running" \
                           ["STOPPED"]="running" \
                           ["SUSPENDED"]="failed" \
                           ["TIMEOUT"]="failed")

echo ${code_to_state[$job_state_code]:-failed}
