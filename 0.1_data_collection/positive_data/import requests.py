import requests
import pandas as pd
import re


#we set the parameters of the query 
url = "https://rest.uniprot.org/uniprotkb/stream"
params = {
    "query": "(reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) "
             "AND (length:[40 TO *]) AND (existence:1) AND (ft_signal_exp:*)",
    "format": "tsv",
    "fields": "accession,organism_name,lineage,length,ft_signal,sequence",
}

#get method from the requests library to save the results of the query. We also use raise_for_status to avoid infinite loop in case of issues
r = requests.get(url, params = params)
r.raise_for_status()

#Create a ew tsv file with the requested parameters
with open("uniprot_results.tsv", "w") as f:
    f.write(r.text)

df = pd.read_csv("uniprot_results.tsv", sep = "\t")


#We use compile method from re library to isolate a specific part of the string which contains the numbers corresponding to the signal peptide length
pat = re.compile(r'^SIGNAL\s+(\d+)\.\.(\d+)\b')

#We create a function to filter the dataset, removing entries with unknown or absent cleavage site
def sp_ok(v):
    if pd.isna(v):
        return False
    t = str(v).strip()
    low = t.lower()
    if "not cleaved" in low or "unknown" in low or any(c in t for c in "?<>"):
        return False
    m = pat.match(t)
    return bool(m) and int(m.group(2)) >= 14

#Apply the function sp_ok at every cell of the column col. Col is a column of True/False, to make us understand if a entry is cleaved (true) or not cleaved or unknown (false)
col = df["Signal peptide"].apply(sp_ok)

#We save the clean columns to the clean variable and we create a new tsv file with all the excluded entries (just to verify if the script worked)
clean = df[col].copy()
df[~col].to_csv("positive_SP_excluded.tsv", sep="\t", index=False)

#for ordered output as request

#We create a fucntion to split the row taxonomic lineage at commas and clean the spaces: the output will be a list like Metazoa (kingdom). Then it rebuilds the set of clean names of the kingdoms
def kingdom(s):
    parts = [p.strip() for p in str(s).split(',')]
    names = {re.sub(r'\s*\(.*\)$', '', p) for p in parts}
    if 'Metazoa' in names:         return 'Metazoa'
    if 'Fungi' in names:           return 'Fungi'
    if 'Viridiplantae' in names:   return 'Plants'
    return 'Other'

#We create a function that reuse the same pattern of the filter to extract the final position of the SP, that is the cleavage site
def cleavage_site(s):
    m = pat.match(str(s).strip())
    return int(m.group(2)) if m else None

#Create a dataframe with exactly the 5 columns requested by the slides, in the correct order
out = pd.DataFrame({
    'accession':     clean['Entry'],
    'organism':      clean['Organism'],
    'kingdom':       clean['Taxonomic lineage'].apply(kingdom),
    'length':        clean['Length'].astype(int),
    'cleavage_site': clean['Signal peptide'].apply(cleavage_site),
    'sequence':      clean['Sequence']     
})

out.to_csv('positive_data.tsv', sep='\t', index=False)

#Build the fasta file iterating one row at a time, creating a tuple of 2 elements: one is the sequence name and the other is the sequence itself
with open("positive_SP_clean.fasta", "w") as f:
    for _, row in clean.iterrows():
        f.write(f">{row['Entry']}\n{row['Sequence']}\n")
