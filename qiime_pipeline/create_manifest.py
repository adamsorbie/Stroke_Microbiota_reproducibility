import pandas as pd 
import numpy as np
import sys

manifest = pd.read_csv(sys.argv[1], names=["sample-id", "absolute-filepath"])

manifest['direction'] = np.where(manifest["sample-id"].str.contains("R1"), 'forward', 'reverse')

manifest.to_csv(sys.argv[2],index=False)

