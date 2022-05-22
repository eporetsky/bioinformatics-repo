
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('text', action='store', type=str, help='PRJNA accession number')
args = parser.parse_args()

# Based on the EBI ENA tutorial on accessing their data
# https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html
# Change Study Accession (PRJNAXXXXXX) to get the 
# Change "fields=all" to get all columns in the metadata table
res = "read_run"
fie = "study_accession,tax_id,scientific_name,instrument_model,library_strategy,read_count,run_alias,sample_alias,fastq_ftp,fastq_md5"

# Generate the url that will containt the study accession information once opened
url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession={0}&result={1}&fields={2}".format(args.text, res, fie)
# Open the URL with pandas and load it as a dataframe
df = pd.read_csv(url, sep="\t")
df.to_csv(args.text+".csv", index=False)