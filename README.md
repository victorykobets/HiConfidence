# HiConfidence
This package helps to construct confidence matrices for Hi-C replicates

## Usage
1. Generate Confidence matrix for the pair of replicates:
```
calculate_confidence(Hi-C_matrix_rep1,Hi-C_matrix_rep2,k_parameter)
```
2. Choose the best parameter **k** based on correlation of IS with *window_size*:
```
choose_k(Hi-C_reps_file_names)
```
3. Construct chromatin profile in the area around averaged TADs boundary using correction for replicates' confidence:
```
average_tad(Hi-C_reps_file_names, k_parameter, tads_positions, area_size)
```
4. Generate Hi-C maps in lower resolution using correction for replicates' confidence:
```
downsampling(Hi-C_reps_file_names, k_parameter, factor, output_file_name)
```

See test_example.py for more examples.
