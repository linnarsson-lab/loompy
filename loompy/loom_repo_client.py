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
	def __init__(self, endpoint: str) -> None:
		self.endpoint = endpoint
	
	def search(self, group: str) -> Dict[str, LoomRepo]:
		return {
			"repos": [
				("slinnarsson", "Lab of Sten Linnarsson at Karolinska Institute"),
				("amitlab", "Ido Amit, Weizmann"),
				("aerts", "Stein Aerts"),
				("niryosef", "Nir Yosef at Berkeley")
			],
			""
		}
