#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fixture_dir="${repo_root}/tests/fixtures/gpu_small"

mlipper_binary="${MLIPPER_TEST_BINARY:-${repo_root}/MLIPPER}"
gpu_id="${MLIPPER_TEST_GPU_ID:-0}"
ref_msa="${MLIPPER_TEST_REF_MSA:-${fixture_dir}/reference.fa}"
query_msa="${MLIPPER_TEST_QUERY_MSA:-${fixture_dir}/query.fa}"
backbone_tree="${MLIPPER_TEST_BACKBONE_TREE:-${fixture_dir}/backbone.nwk}"
model="${MLIPPER_TEST_MODEL:-${fixture_dir}/model.bestModel}"
dnc_msa="${MLIPPER_TEST_DNC_MSA:-${fixture_dir}/full.fa}"
expected_jplace="${MLIPPER_TEST_EXPECTED_JPLACE:-${fixture_dir}/expected.jplace}"
expected_dnc_tree="${MLIPPER_TEST_EXPECTED_DNC_TREE:-${fixture_dir}/expected.dnc.nwk}"
test_dir="$(mktemp -d)"
trap 'rm -rf -- "${test_dir}"' EXIT

"${mlipper_binary}" \
    --tree-alignment "${ref_msa}" \
    --query-alignment "${query_msa}" \
    --tree "${backbone_tree}" \
    --best-model "${model}" \
    --jplace-out "${test_dir}/actual.jplace" \
    --gpu-id "${gpu_id}" \
    > "${test_dir}/placement.log" 2>&1

"${mlipper_binary}" \
    --divide-and-conquer \
    --tree-alignment "${dnc_msa}" \
    --best-model "${model}" \
    --dipper-starting-tree-mode "${MLIPPER_TEST_DIPPER_TREE_MODE:-nj-placement}" \
    --divide-and-conquer-core-edges "${MLIPPER_TEST_DNC_CORE_EDGES:-2}" \
    --divide-and-conquer-tip-budget "${MLIPPER_TEST_DNC_TIP_BUDGET:-6}" \
    --divide-and-conquer-sweeps "${MLIPPER_TEST_DNC_SWEEPS:-1}" \
    --gpu-id "${gpu_id}" \
    --write-tree "${test_dir}/actual.dnc.nwk" \
    > "${test_dir}/dnc.log" 2>&1

python3 "${repo_root}/tests/compare_workflow_outputs.py" \
    --expected-jplace "${expected_jplace}" \
    --actual-jplace "${test_dir}/actual.jplace" \
    --expected-tree "${expected_dnc_tree}" \
    --actual-tree "${test_dir}/actual.dnc.nwk"

grep -q "NNI likelihood audit:" "${test_dir}/dnc.log"
grep -q "MLIPPER D&C final optimization done:" "${test_dir}/dnc.log"
echo "workflow_regression: PASS"
