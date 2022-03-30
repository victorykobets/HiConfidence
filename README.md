# ReSiE
This package helps to construct confidence matrices for Hi-C replicates

## Usage
1. Generate Confidence matrix for the pair of replicates:
```
get_confidence(matrix_1,matrix_2,k)
```
2. Choose the best parameter **k** based on correlation of IS with *window_size*:
```
choose_k(file_names=[*replicates], max_k = 5, window_size)
```
3. Construct chromatin profile in the area around averaged TADs boundary using correction for replicates' confidence:
```
average_tad(file_names=[*replicates], k, tads_positions=[chr,start,end...], area_size)
```
4. Generate Hi-C maps in lower resolution using correction for replicates' confidence:
```
downsampling(file_names=[*replicates], k=2, multiplier=2, assembly_version, generate_files=[*new_files_names])
```
