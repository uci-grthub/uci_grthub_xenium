#!/usr/bin/env python
import pandas as pd
import scanpy as sc
from pathlib import Path

adata = sc.read(Path("output/scanpy/all_harmony_integated.h5ad")) 
                         
aggregated = sc.get.aggregate(adata, by="batch", func=["mean", "count_nonzero"])
