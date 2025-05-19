# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns
import glob as glob
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator, ScalarFormatter, MultipleLocator
import sys

# Check arguments
if len(sys.argv) < 3:
    print("Usage: python3 processing.py <input_file1> <input_file2> ... <output_file>")
    sys.exit(1)

# Define input and output file 
input = sys.argv[1:-1]
output = sys.argv[-1]

# Initialize an empty DataFrame
df = pd.DataFrame()

#### PROCESSING FILES ####
# Loop through all input files and process them
for i in input:
    print(f"Processing {i}")
    # Check if the file has a .csv extension; if not, exit with an error message
    if not i.endswith(".csv"):
        print("Input files must be in .csv format")
        sys.exit(1)
    else:
        # Read the CSV file into a DataFrame with custom column names and data types
        dataset = pd.read_csv(i, sep=",", names=["Platform", "Length", "Quality"],
                              dtype={"Platform": str, "Length": float, "Quality": float})
        
        # Group data by Platform, Length, and Quality, and count occurrences
        agg = dataset.groupby(['Platform', 'Length', 'Quality']).size().reset_index(name='Counts')
        
        # Save the aggregated data to the output file
        agg.to_csv(output, index=False)
    
    print(f"Concatenating {i}")
    # Concatenate the current aggregated data 
    df = pd.concat([df, agg])

df.to_csv(output, index=False)

# Read the aggregated data from the output file
df = pd.read_csv(output, sep=",", header=0)


### LENGTH ASSESSMENT ###
# Filter the data to include only ONT and Hifi platforms
lr = df[(df['Platform'] == 'ONT') | (df['Platform'] == 'Hifi')]

# Apply a log10 transformation to the 'Length' column for better visualization
lr['Length'] = np.log10(lr['Length'])

# Create a new figure with specific size
plt.figure(figsize=(9, 8))

# Define color mappings for different platforms
colors = {'ONT': '#118AB2', 'Hifi': '#EF476F', 'illumina': '#154CAC', 'AVITI': '#06D6A0', '10x': '#FFD166'}

# Plot a histogram of 'Length' for the filtered data, with customization
ax = sns.histplot(data=lr, x='Length', hue='Platform', edgecolor='white',
                  hue_order=['Hifi','ONT'],weights='Counts', fill=True, palette=colors, bins=100,
                  legend=False, alpha=1)

# Highlight a region in the plot (for example, between x-values 2 and 2.55)
ax.axvspan(2, 2.55, color='#154CAC', alpha=1)

# Set the x and y axis limits and labels
plt.xlim(1.8, 5)
plt.yticks(ticks=[50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000, 500000],
           labels=['50', '100', '150', '200', '250', '300', '350', '400', '450', '500'], fontsize=26)
plt.xticks(ticks=[2, 3, 4, 5, 6], labels=[100, 1000, "10K", "100K", "1000K"], fontsize=26)
plt.xlabel('Read Length (bp)', fontsize=26)
plt.ylabel('Number of Reads (K)', fontsize=26)
plt.rcParams['svg.fonttype'] = 'none'
plt.savefig('Figure_2b.svg', bbox_inches='tight')
