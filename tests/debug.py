
import loompy
import numpy as np

global_attrs = {
	"SampleID": "kjhlkjh",
	"AmbientUMIs": 2983470867,
	"MinimumCellUMIs": 2343523,
	"Saturation": 0.5903294892
}
loompy.create("/Users/stelin/kallisto_GRCh38/10X_17_029b.loom", np.zeros((100, 100)), {}, {}, file_attrs=global_attrs)
