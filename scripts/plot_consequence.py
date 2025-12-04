#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def plot_consequence(df: pd.DataFrame, output_path: str) -> None:
    consequence_counts = df["Consequence"].value_counts()

    plt.figure(figsize=(10, 6))
    consequence_counts.plot(kind='bar')
    plt.title('Variant Consequence Distribution')
    plt.xlabel('Consequence Type')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def main():
    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <input_tsv> <output_png>\n")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    df = pd.read_table(input_path, sep='\t')
    plot_consequence(df, output_path)

if __name__ == '__main__':
    main()