# %%
from dask.distributed import Client
client = Client('tcp://127.0.0.1:8786')
client

# %%
import intake
intake.output_notebook()  # enables source.plot API

import hvplot.xarray

# %% URL to github repo is NOT working because it's PRIVATE; use local file.
#icecaps = intake.open_catalog('https://raw.githubusercontent.com/vonw/icecaps-data/main/intake/icecaps.yaml?token=GHSAT0AAAAAACELOP5T3ECL33J4UF5AOPPQZNIFI5Q')
icecaps = intake.open_catalog('/Users/vonw/work/software/icecaps-data/intake/icecaps.yaml')
list(icecaps)

# %%
sfcmet = icecaps.sfcmet.to_dask()
sfcmet

# %% This may take a while, depending on internet access to U. Leeds server and local computing power...
T2m = sfcmet['temperature at 2 meters'].resample(time='30Min').mean().compute()

# %%
