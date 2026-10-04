#!/usr/bin/env bash
# sync_skills.sh
# Synchronizes agent skills between the git repository (skills/) and the user configuration (~/.gemini/config/skills/).
#
# Usage:
#   ./scripts/sync_skills.sh [--check | --to-system | --to-repo]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_SKILLS="${REPO_ROOT}/skills"
USER_SKILLS="${HOME}/.gemini/config/skills"

MODE="${1:---to-system}"

echo "=================================================================="
echo "SKILL SYNCHRONIZATION UTILITY"
echo "  Repo Skills Path : ${REPO_SKILLS}"
echo "  User Skills Path : ${USER_SKILLS}"
echo "  Sync Mode        : ${MODE}"
echo "=================================================================="

mkdir -p "${USER_SKILLS}"
mkdir -p "${REPO_SKILLS}"

if [ "${MODE}" == "--check" ]; then
    echo "[Drift Check] Comparing repository skills with user config..."
    DRIFT_FOUND=0
    for skill_dir in "${REPO_SKILLS}"/*; do
        if [ -d "${skill_dir}" ]; then
            skill_name="$(basename "${skill_dir}")"
            user_skill_dir="${USER_SKILLS}/${skill_name}"
            if [ ! -d "${user_skill_dir}" ]; then
                echo "  [DRIFT] Skill '${skill_name}' exists in repo but missing in user config."
                DRIFT_FOUND=1
            elif ! diff -r -q "${skill_dir}" "${user_skill_dir}" >/dev/null 2>&1; then
                echo "  [DRIFT] Skill '${skill_name}' content differs between repo and user config."
                diff -u "${skill_dir}/SKILL.md" "${user_skill_dir}/SKILL.md" || true
                DRIFT_FOUND=1
            fi
        fi
    done
    if [ ${DRIFT_FOUND} -eq 0 ]; then
        echo "[Drift Check] Zero drift detected. All skills are in perfect synchronization."
    else
        echo "[Drift Check Warning] Drift detected. Run './scripts/sync_skills.sh --to-system' to synchronize."
        exit 1
    fi
elif [ "${MODE}" == "--to-system" ]; then
    echo "[Sync] Copying skills from repository -> user config..."
    for skill_dir in "${REPO_SKILLS}"/*; do
        if [ -d "${skill_dir}" ]; then
            skill_name="$(basename "${skill_dir}")"
            target_dir="${USER_SKILLS}/${skill_name}"
            mkdir -p "${target_dir}"
            cp -r "${skill_dir}"/* "${target_dir}/"
            echo "  [Synced] ${skill_name} -> ${target_dir}/"
        fi
    done
    echo "[Sync Complete] User config skills are up to date."
elif [ "${MODE}" == "--to-repo" ]; then
    echo "[Sync] Copying skills from user config -> repository..."
    for skill_dir in "${USER_SKILLS}"/*; do
        if [ -d "${skill_dir}" ]; then
            skill_name="$(basename "${skill_dir}")"
            target_dir="${REPO_SKILLS}/${skill_name}"
            mkdir -p "${target_dir}"
            cp -r "${skill_dir}"/* "${target_dir}/"
            echo "  [Synced] ${skill_name} -> ${target_dir}/"
        fi
    done
    echo "[Sync Complete] Repository skills are up to date."
else
    echo "Unknown option: ${MODE}. Valid options: --check, --to-system, --to-repo"
    exit 1
fi
