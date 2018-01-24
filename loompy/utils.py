import logging
from typing import *
from inspect import currentframe, getouterframes
import datetime


def deprecated(message: str) -> None:
	frameinfo = getouterframes(currentframe())
	logging.warn(f"╭── " + message)
	logging.warn(f"╰──> at {frameinfo[2].filename}, line {frameinfo[2].lineno}")


def timestamp() -> str:
	return datetime.datetime.utcnow().strftime("%Y%m%dT%H%M%S.%fZ")
