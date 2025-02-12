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
        print("Input files must be in .tsv format")
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

# Read the aggregated data from the output file
df = pd.read_csv(output, sep=",", header=0)


### QUALITY ASSESSMENT ####
# Filter the DataFrame to include only rows where 'Quality' is less than 50 to  remove reads of high quality due to  short length
agg2 = df[df['Quality'] < 50]

# Define color mappings for different platforms
colors = {'ONT': '#118AB2', 'Hifi': '#EF476F', 'illumina': '#154CAC', 'AVITI': '#06D6A0', '10x': '#FFD166'}
# Define the order of platforms for plotting
order = ["ONT", "Hifi", "AVITI", "10x", "illumina"]

print("Start plotting")
# Create a FacetGrid for plotting histograms for each platform, with customization
g = sns.FacetGrid(agg2, row="Platform", aspect=4, height=1, sharex=True, sharey=False, row_order=order)
g.map_dataframe(sns.histplot, x="Quality", weights="Counts", binwidth=1, hue="Platform", palette=colors, edgecolor='white', legend=True, alpha=1)
# Remove titles for each row (platform)
g.set_titles(row_template="")
# Customize the y-axis labels for each platform
for ax, platform in zip(g.axes.flat, order):
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))  # Set the number of bins on the y-axis
    if platform == 'ONT':
        ax.set_yticklabels(["0", "150", "300"])
        ax.set_ylabel('')
    elif platform == 'Hifi':
        ax.set_yticklabels(["0", "30", "60"])
        ax.set_ylabel('')
    elif platform == 'illumina':
        ax.set_yticklabels(["0", "2000", "4000"])
        ax.set_ylabel('')
    elif platform == '10x':
        ax.set_yticklabels(["0", "15000", "30000"])
        ax.set_ylabel('')
    else:
        ax.set_yticklabels(["0", "600000", "120000"])
        ax.set_ylabel('')
g.figure.text(0.1, 0.5, 'Number of Reads (K)', va='center', rotation=90, fontsize=12)
g.set_xlabels('Quality Score')
g.fig.subplots_adjust(hspace=0.1)
g.fig.set_size_inches(6, 6)
plt.rcParams['svg.fonttype'] = 'none'
plt.savefig('Quality.svg', bbox_inches='tight')

