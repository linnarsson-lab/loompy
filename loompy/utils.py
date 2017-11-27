import logging
from typing import *
from inspect import currentframe, getouterframes


def deprecated(message: str) -> None:
	frameinfo = getouterframes(currentframe())
	logging.warn(f"╭── " + message)
	logging.warn(f"╰──> at {frameinfo[2].filename}, line {frameinfo[2].lineno}")
