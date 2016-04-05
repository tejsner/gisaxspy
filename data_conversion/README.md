# Converting files from SPEC
To convert files from a SPEC experiement, use `convert_from_spec.py` as follows:

```
from convert_from_spec import *
convert_from_spec('pspb05_all', 'ex_gisaxs.csv', 31, 547, samx_col=4)
```

This takes `pspb05_all` as input, `ex_gisaxs.csv` as output and extracts data from image 32 to 547. When doing this conversion it is important that there are no "breaks" in the experiment (XR, lineups, etc) as the parser is rather simple. The `samx_col` parameter specifies what column (in the line beginning with `#P1`) that contains information about the x-axis of the experiment.
