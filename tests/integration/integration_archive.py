#!/usr/bin/env python
import argparse
import shutil
from os import path
import os

import logging
logger = logging.getLogger("integration_archive")
logging.basicConfig(level=logging.INFO)

def save_results(archive, results, branch, revision, force=False):
    revision_path = path.join(archive, "revisions", revision)

    if path.lexists(revision_path):
        if force:
            logger.info("Removing existing result: %s", revision_path)
            shutil.rmtree(revision_path)
        else:
            raise ValueError("Archive results for revision %s already exists: %s" % (revision, revision_path))

    logger.info("Saving result at revision path: %s", revision_path)
    shutil.copytree(results, revision_path)

    branch_path = path.join(archive, branch)
    if path.lexists(branch_path):
        logger.info("Removing exsiting branch pointer: %s targeting: %s", branch_path, os.readlink(branch_path))
        os.remove(branch_path)

    logger.info("Saving branch pointer: %s targeting: %s", branch_path, path.join("revision", revision))
    os.symlink(path.join("revisions", revision), branch_path)

def save_results_main(args):
    save_results(args.archive, args.results, args.branch, args.revision, force=args.force)

parser = argparse.ArgumentParser(prog="integration_archive.py")

subparsers = parser.add_subparsers()

save_parser = subparsers.add_parser("save", help="Save integration test results in archive")
save_parser.add_argument("--force", default=False, action='store_true', help="Overwrite results if revision is already archived.")
save_parser.add_argument("archive", help="Integration archive directory.")
save_parser.add_argument("results", help="Result directory.")
save_parser.add_argument("branch", help="Result source branch.")
save_parser.add_argument("revision", help="Result source revision.")
save_parser.set_defaults(func=save_results_main)

args = parser.parse_args()
args.func(args)
