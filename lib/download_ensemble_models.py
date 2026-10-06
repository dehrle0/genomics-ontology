#!/usr/bin/env python3
"""
download_ensemble_models.py
Downloads the specialized ensemble models for the clinical multi-track pipeline:
1. QwQ-32B (Adversarial Review & Conflict Adjudication)
2. Mistral-Small-3.2-24B (Publication Report Writer)
3. Baichuan-M2-32B (Pharmacogenomics & Clinical Drug Response)
4. Bio-Medical-Llama-3.1-8B (Genomic Co-Factor & Literature Validation)
"""

import os
import sys
import time
from huggingface_hub import hf_hub_download

MODELS = [
    {
        "name": "Bio-Medical-Llama-3.1-8B",
        "repo": "mradermacher/Bio-Medical-Llama-3.1-8B-GGUF",
        "filename": "Bio-Medical-Llama-3.1-8B.Q4_K_M.gguf",
        "target_name": "Bio-Medical-Llama-3.1-8B-Q4_K_M.gguf"
    },
    {
        "name": "Mistral-Small-3.2-24B",
        "repo": "unsloth/Mistral-Small-3.2-24B-Instruct-2506-GGUF",
        "filename": "Mistral-Small-3.2-24B-Instruct-2506-Q4_K_M.gguf",
        "target_name": "Mistral-Small-3.2-24B-Instruct-2506-Q4_K_M.gguf"
    },
    {
        "name": "QwQ-32B",
        "repo": "Qwen/QwQ-32B-GGUF",
        "filename": "qwq-32b-q4_k_m.gguf",
        "target_name": "qwq-32b-q4_k_m.gguf"
    },
    {
        "name": "Baichuan-M2-32B",
        "repo": "mradermacher/Baichuan-M2-32B-GGUF",
        "filename": "Baichuan-M2-32B.Q4_K_M.gguf",
        "target_name": "Baichuan-M2-32B-Q4_K_M.gguf"
    }
]

TARGET_DIR = os.path.expanduser("~/ai-infrastructure/llama.cpp/models")
LMSTUDIO_DIR = os.path.expanduser("~/.lmstudio/models/ensemble")

def main():
    os.makedirs(TARGET_DIR, exist_ok=True)
    os.makedirs(LMSTUDIO_DIR, exist_ok=True)

    print("==================================================================")
    print("STARTING CLINICAL ENSEMBLE GGUF MODEL DOWNLOAD SEQUENCE")
    print(f"Target Directory: {TARGET_DIR}")
    print("==================================================================")

    for idx, m in enumerate(MODELS, start=1):
        target_path = os.path.join(TARGET_DIR, m["target_name"])
        symlink_path = os.path.join(LMSTUDIO_DIR, m["target_name"])

        print(f"\n[{idx}/{len(MODELS)}] Checking {m['name']} ({m['repo']})...")
        if os.path.exists(target_path):
            size_gb = os.path.getsize(target_path) / (1024**3)
            print(f"  -> Already exists: {target_path} ({size_gb:.2f} GB). Skipping download.")
        else:
            t0 = time.time()
            print(f"  -> Downloading {m['filename']} from {m['repo']}...")
            try:
                downloaded_file = hf_hub_download(
                    repo_id=m["repo"],
                    filename=m["filename"],
                    local_dir=TARGET_DIR,
                    local_dir_use_symlinks=False
                )
                if m["filename"] != m["target_name"]:
                    actual_downloaded = os.path.join(TARGET_DIR, m["filename"])
                    if os.path.exists(actual_downloaded):
                        os.rename(actual_downloaded, target_path)
                dt = time.time() - t0
                size_gb = os.path.getsize(target_path) / (1024**3)
                print(f"  -> Successfully downloaded {m['name']} in {dt:.1f}s ({size_gb:.2f} GB).")
            except Exception as e:
                print(f"  -> ERROR downloading {m['name']}: {e}")
                sys.exit(1)

        # Ensure symlink in LM Studio directory for auto-discovery
        if not os.path.exists(symlink_path):
            try:
                os.symlink(target_path, symlink_path)
                print(f"  -> Linked to LM Studio: {symlink_path}")
            except Exception as e:
                print(f"  -> Link warning: {e}")

    print("\n==================================================================")
    print("ALL ENSEMBLE MODELS VERIFIED AND READY FOR SEQUENTIAL INFERENCE")
    print("==================================================================")

if __name__ == "__main__":
    main()
