from typing import *
import json

"""
REST API

GET /?query={query}				Return JSON list(s) of organizations, repos and datasets that match query
GET /{org}						Return JSON list of repos for this organization
GET /{org}/{repo}				Return JSON list of datasets for this repo
GET /{org}/{repo}/{dataset}		Return the specific dataset
PUT /{org}/{repo}/{dataset}		Upload a dataset
DELETE /{org}/{repo}/{dataset}	Remove a dataset

# TODO: create, update and delete organizations
#       create, update and delete users
#       assign user to organization (and revoke user from organization)
#       create, update and delete repository


PYTHON API

>>> client = LoomRepoClient("http://loompy.org/api")
>>> client.search("linnarsson") # Default is to return a pretty string
linnarsson-lab (organization)		Lab of Sten Linnarsson at Karolinska Institute

>>> client["linnarsson-lab"].repos()		# List of all repos under this organization
linnarsson-lab (organization)		Lab of Sten Linnarsson at Karolinska Institute
-----------------------------------------------------------------------------------
cortex-2015							Data for Zeisel et al., Science 2015
mouse-brain-2018					Data for Zeisel et al., Cell 2018
sympathetic							Data for Furlan et al., Nature Neuroscience 2017 

>>> client["linnarsson-lab"]["cortex-2015"]		# List all datasets under this repo
linnarsson-lab/cortex-2015 (repository)		Data for Zeisel et al., Science 2015
-----------------------------------------------------------------------------------
cortex-2015.loom							27995 genes and 3005 cells

>>> client["linnarsson-lab"]["cortex-2015"]["cortex-2015.loom"]  # Show info for this dataset
linnarsson-lab/cortex-2015/cortex-2015.loom (file)		Data for Zeisel et al., Science 2015
-----------------------------------------------------------------------------------
Description:	ljlkj hlkjhlkhklh 
      Shape:	(27995, 3005)
	Journal:	Science
	   Year:	2015

>>> client["linnarsson-lab"]["cortex-2015"]["cortex-2015.loom"].download()
(shows a progress bar while downloading)

>>> client["linnarsson-lab/cortex-2015/cortex-2015.loom"].download()  # Equivalent syntax

BASH API

$ loom search linnarsson
linnarsson-lab (organization)		Lab of Sten Linnarsson at Karolinska Institute

$ loom list linnarsson-lab
linnarsson-lab (organization)		Lab of Sten Linnarsson at Karolinska Institute
-----------------------------------------------------------------------------------
cortex-2015							Data for Zeisel et al., Science 2015
mouse-brain-2018					Data for Zeisel et al., Cell 2018
sympathetic							Data for Furlan et al., Nature Neuroscience 2017 

$ loom list "linnarsson-lab/cortex-2015"
linnarsson-lab/cortex-2015 (repository)		Data for Zeisel et al., Science 2015
-----------------------------------------------------------------------------------
cortex-2015.loom							27995 genes and 3005 cells

$ loom show linnarsson-lab/cortex-2015/cortex-2015.loom
linnarsson-lab/cortex-2015/cortex-2015.loom (file)		Data for Zeisel et al., Science 2015
-----------------------------------------------------------------------------------
Description:	ljlkj hlkjhlkhklh 
      Shape:	(27995, 3005)
	Journal:	Science
	   Year:	2015

$ loom get "linnarsson-lab/cortex-2015/cortex-2015.loom"
(shows a progress bar while downloading)

"""


class SearchResult:
	def __init__(self, org: str, org_desc: str, repo: str = None, repo_desc: str = None, dataset: str = None, dataset_desc: str = None) -> None:
		self.org = org
		self.org_desc = org_desc
		self.repo = repo
		self.repo_desc = repo_desc
		self.dataset = dataset
		self.dataset_desc = dataset_desc

class LoomRepoClient:
	def __init__(self, endpoint: str) -> None:
		self.endpoint = endpoint
		self.organizations = dict()

	def search(self, group: str) -> Dict[str, List[Tuple]]:
		return [
			SearchResult("linnarsson-lab", "Lab of Sten Linnarsson at Karolinska Institute"),
			SearchResult("nir-yosef-lab", "Lab of Nir Yosef at Berkeley")
			SearchResult("linnarsson-lab", "Lab of Sten Linnarsson at Karolinska Institute", "mouse-brain", "Mouse Brain Architecture Cell 2018", "L5_All.loom", "All cells that passed QC")
		]
