from typing import *


def to_html(ds: Any) -> str:
	"""
	Return an HTML representation of the loom file or view, showing the upper-left 10x10 corner.
	"""
	rm = min(10, ds.shape[0])
	cm = min(10, ds.shape[1])
	html = "<p>"
	if ds.attrs.__contains__("title"):
		html += "<strong>" + ds.attrs["title"] + "</strong> "
	html += f"{ds.shape[0]} rows, {ds.shape[1]} columns, {len(ds.layers)} layer{'s' if len(ds.layers) > 1 else ''}<br/>(showing up to 10x10)<br/>"
	html += ds.filename + "<br/>"
	for (name, val) in ds.attrs.items():
		html += f"name: <em>{val}</em><br/>"
	html += "<table>"
	# Emit column attributes
	for ca in ds.col_attrs.keys():
		html += "<tr>"
		for ra in ds.row_attrs.keys():
			html += "<td>&nbsp;</td>"  # Space for row attrs
		html += "<td><strong>" + ca + "</strong></td>"  # Col attr name
		for v in ds.col_attrs[ca][:cm]:
			html += "<td>" + str(v) + "</td>"
		if ds.shape[1] > cm:
			html += "<td>...</td>"
		html += "</tr>"

	# Emit row attribute names
	html += "<tr>"
	for ra in ds.row_attrs.keys():
		html += "<td><strong>" + ra + "</strong></td>"  # Row attr name
	html += "<td>&nbsp;</td>"  # Space for col attrs
	for v in range(cm):
		html += "<td>&nbsp;</td>"
	if ds.shape[1] > cm:
		html += "<td>...</td>"
	html += "</tr>"

	# Emit row attr values and matrix values
	for row in range(rm):
		html += "<tr>"
		for ra in ds.row_attrs.keys():
			html += "<td>" + str(ds.row_attrs[ra][row]) + "</td>"
		html += "<td>&nbsp;</td>"  # Space for col attrs

		for v in ds[row, :cm]:
			html += "<td>" + str(v) + "</td>"
		if ds.shape[1] > cm:
			html += "<td>...</td>"
		html += "</tr>"
	# Emit ellipses
	if ds.shape[0] > rm:
		html += "<tr>"
		for v in range(rm + 1 + len(ds.row_attrs.keys())):
			html += "<td>...</td>"
		if ds.shape[1] > cm:
			html += "<td>...</td>"
		html += "</tr>"
	html += "</table>"
	return html
