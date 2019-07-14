import loompy
import os
loompy.combine_faster([os.path.join("/Users/stelin/loom/", x) for x in ("10X160_1.loom", "10X161_1.loom")], "/Users/stelin/combined.loom")