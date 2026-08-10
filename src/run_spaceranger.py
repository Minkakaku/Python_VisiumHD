#!/usr/bin/env python
"""从 YAML 配置构建并运行 Space Ranger count 命令。"""

import argparse
import shlex
import subprocess
from pathlib import Path

import yaml


# 功能：把 YAML 项目参数转换成不含硬编码路径的 spaceranger count 参数列表。
def build_spaceranger_command(config):
    command = ["spaceranger", "count", f"--id={config['id']}"]
    scalar_options = {
        "transcriptome": "transcriptome",
        "sample": "sample",
        "image": "image",
        "slide": "slide",
        "area": "area",
        "probe_set": "probe-set",
        "cytaimage": "cytaimage",
        "loupe_alignment": "loupe-alignment",
        "localcores": "localcores",
        "localmem": "localmem",
    }
    for key, option in scalar_options.items():
        if config.get(key) not in (None, ""):
            command.append(f"--{option}={config[key]}")
    fastqs = config.get("fastqs", [])
    if fastqs:
        command.append(f"--fastqs={','.join(map(str, fastqs))}")
    if "create_bam" in config:
        command.append(f"--create-bam={str(bool(config['create_bam'])).lower()}")
    return command


# 功能：运行 Space Ranger 并把完整命令与标准输出记录到日志。
def run_spaceranger(config_file, log_file=None, dry_run=False):
    config_file = Path(config_file)
    config = yaml.safe_load(config_file.read_text(encoding="utf-8"))
    command = build_spaceranger_command(config)
    printable = shlex.join(command)
    if dry_run:
        print(printable)
        return command
    log_file = Path(log_file or f"logs/{config['id']}.spaceranger.log")
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("w", encoding="utf-8") as handle:
        handle.write(f"COMMAND\n{printable}\n\nOUTPUT\n")
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT, text=True)
    return command


# 功能：解析 Space Ranger 独立入口参数。
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML configuration file")
    parser.add_argument("--log", help="log file; default logs/<id>.spaceranger.log")
    parser.add_argument("--dry-run", action="store_true", help="print command without running")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_spaceranger(args.config, args.log, args.dry_run)
