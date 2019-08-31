import logging
import os
import sys
from typing import List

import click

from ._version import __version__
from .bus_file import create_from_fastq


@click.group()
@click.option('--show-message/--hide-message', default=True)
@click.option('--verbosity', default="info", type=click.Choice(['error', 'warning', 'info', 'debug']))
def cli(show_message: bool = True, verbosity: str = "info") -> None:
	level = {"error": 40, "warning": 30, "info": 20, "debug": 10}[verbosity]
	logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=level)
	logging.captureWarnings(True)

	if show_message:
		print(f"Loompy v{__version__} by Linnarsson Lab 🌸 (http://linnarssonlab.org & http://loompy.org)")
		print()


@cli.command()
@click.argument('loomfile', required=True)
@click.argument('sampleid', required=True)
@click.argument('indexdir', required=True)
@click.argument('metadatafile', required=True)
@click.argument('fastqs', required=True, nargs=-1)
@click.option('--threads', default=os.cpu_count(), help="Number of threads to use")
def fromfq(loomfile: str, sampleid: str, indexdir: str, metadatafile: str, threads: int, fastqs: List[str]) -> None:
	logging.info(f"Using {threads} threads.")
	create_from_fastq(loomfile, sampleid, list(fastqs), indexdir, metadatafile, threads)
