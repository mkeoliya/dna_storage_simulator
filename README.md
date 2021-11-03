# DNA Storage Simulator

This monorepository contains code for a research project by Mayank Keoliya and supervised by Djordje Jevdjic, that analyzes and simulates DNA storage. 

The structure is as follows:

- `CodeReconstruction`: Forked from [CodeReconstruction](https://github.com/omersabary/CodeReconstruction), with modifications to aid testing.
    - `real_data_clustered.txt`: Stores real data in the form:
        ```
        [original strand][\n]
        *****************************[\n]
        [copy][\n]
        [copy][\n]
        ...
        [copy][\n]
        [\n]
        [\n]
        [original strand][\n]
        *****************************[\n]
        [copy][\n]
        [copy][\n]
        ...
        [copy][\n]
        ...
        ```
    - `synth_data_clustered.txt`: Stores synthetic data in the same form as `real_data_clustered`, can be generated via the `noisy.py` module, or [DNASimulator](https://github.com/gadihh/DNASimulator)
    - `compare.sh`: Bash script that runs reconstruction algorithms on `real_data_clustered` and `synth_data_clustered` 

- `Scripts`: Contains utility scripts
    - `get_ground_from_clustered.py`: Parses files of the same form as `real_data_clustered.txt` above, and generates a file `strands.txt` containing the original strands only. Run using `python get_ground_from_clustered.py`
    - `noisy.py`: Naive simulator that takes in a `strands.txt` file, sequencing coverage, and error probabilities as input and generates noisy copies of multiple clusters in the same form as `real_data_clustered.txt` 