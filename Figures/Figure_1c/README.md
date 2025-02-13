### Generate the Data for Quality and Length Distribution Plot

To generate the datasets for running the **quality and length distribution plots**, follow the instructions below:

1. **Include the technology name** in the file name, separated by **dots** before the file extension.  

   **Example:**  
   `Dataset_101.HiFi.fastq`

2. **Run the `rdeval -qa` flag** for all read files to perform the analysis:

    ```bash
    for file in *.[fasta|fastq|fasta.gz|fastq.gz|bam|cram]; do
        rdeval "$file" -qa > "${file}.ql"
    done
    ```

3. **Generate the technology column** in each dataset by running the following script:

    ```bash
    for i in *.ql; do
        # Extract the penultimate field (technology) from the filename
        name=$(awk -F. '{print $(NF-1)}' <<< "$i")
        
        # Create CSV files with technology and data columns
        awk -v OFS="," -v n=$name '{print n,$1,$2}' "$i" >"${i}.csv"
    done
    ```

    This will generate CSV files with the required format to run the `quality_assessment.py` and `length_assessment.py` scripts. 

4. **In the directory containing the CSV files**, execute the following command to call the required Python script:

    ```bash
    python3 [quality|length].py <input_file1> <input_file2> ... <output_file>
    ```

### Output ###
 **it will generate the plot in svg format and  concatanete dataset.csv format**  
