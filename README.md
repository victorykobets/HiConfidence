# ReSiE
This package helps to construct confidence matrices for Hi-C replicates

## Usage
1. Construct chromatin profile in the area around averaged TADs boundary using correction for replicates confidence:
```
average\_tad(file\_names=\[*replicates], k=2, tads\_positions=\[chr,start,end], area_size=10)
```
2. Generate Hi-C maps in lower resolution using correction for replicates confidence:
```
downsampling(file\_names=\[\*replicates], k=2, multiplier=2, assembly\_v, generate\_files=\[*new\_file\_names])
```
