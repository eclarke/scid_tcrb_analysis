# README
Erik Clarke  

The main document is `report.Rmd`. Supporting code is located in `lib\`.

## How to add new runs to the pipeline:
1. Download the data from Adaptive using the Sample Export. This should be in a .zip format containing a .tsv file for each sample.

2. Download the sample summary data from Adaptive using Sample Summary. This should be in a .tsv format with one row per sample.

3. Append the sample summary data to the existing sample summary .tsv file, after verifying that the columns are the same: 

    ```bash 
    $ tail -n +2 sampleStats.some_timestamp.tsv >> sampleStats.all.tsv
    ```

4. Place the unzipped .tsv files in the `data/tsv` folder.

5. Make sure the `seq_files` points to the `data/tsv` folder, and the `seq_stats_file` points to the complete sample summary file.


6. Add the ng of DNA for each sample to the sampleDNA.xlsx file and export to a tab-separated file (sampleDNA.txt). The ng of DNA can be found in Adaptive's quant reports or from our own in-house measurements before the samples are sent out.
