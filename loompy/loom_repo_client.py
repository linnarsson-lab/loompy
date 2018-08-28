from typing import *
import json

"""
REST API

GET /?query={query}					Return JSON list(s) of groups, repos and datasets that match query
GET /{group}						Return JSON list of repos for this group
GET /{group}/{repo}					Return JSON list of datasets for this repo
GET /{group}/{repo}/{dataset}		Return the specific dataset
PUT /{group}/{repo}/{dataset}		Upload a dataset
DELETE /{group}/{repo}/{dataset}	Remove a dataset
"""


class LoomRepo:
	def __init__(self, group: str, repo: str) -> None:
		self.group = group
		self.repo = repo


class LoomRepoClient:
	def __init__(self, endpoint: str = "http://loom.linnarssonlab.org") -> None:
		self.endpoint = endpoint
	
	def list_repos_for_group(self, group: str) -> Dict[str, LoomRepo]:
		return [LoomRepo(group, )