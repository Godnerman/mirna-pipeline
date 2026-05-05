
import pandas as pd
import os
import re

input_files = snakemake.input
output_file = snakemake.output[0]
samples_csv = snakemake.params.samples_csv

metadata = pd.read_csv(samples_csv)
metadata = metadata.drop_duplicates(subset="gsm")

def parse_mrd(filepath):
    counts = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            m = re.match(r'^(hsa-\S+)\s+read count\s+(\d+)', line)
            if m:
                mirna = m.group(1)
                count = int(m.group(2))
                counts[mirna] = counts.get(mirna, 0) + count
    return counts

all_counts = {}
all_gsms = []

for f in input_files:
    gsm = f.split("/")[-2]
    all_gsms.append(gsm)
    try:
        counts = parse_mrd(f)
        all_counts[gsm] = counts
        print(f"Parsed {gsm}: {len(counts)} miRNAs detected")
    except Exception as e:
        print(f"Warning: error parsing {f}: {e}")
        all_counts[gsm] = {}

all_mirnas = sorted(set(
    mirna for counts in all_counts.values() for mirna in counts
))

print(f"Total unique miRNAs: {len(all_mirnas)}")

matrix = {}
for gsm in all_gsms:
    matrix[gsm] = {m: all_counts[gsm].get(m, 0) for m in all_mirnas}

df = pd.DataFrame(matrix, index=all_mirnas)
df = df.fillna(0).astype(int)
print(f"Matrix: {df.shape[0]} miRNAs x {df.shape[1]} samples")

os.makedirs(os.path.dirname(output_file), exist_ok=True)
df.to_csv(output_file)
print(f"Saved to: {output_file}")
