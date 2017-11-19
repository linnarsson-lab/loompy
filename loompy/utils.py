import logging
from typing import *
from inspect import currentframe, getouterframes


def deprecated(message: str) -> None:
	frameinfo = getouterframes(currentframe())
	logging.warn(f"At {frameinfo[2].filename}, line {frameinfo[2].lineno}:")
	logging.warn(message)
