import pandas as pd
from itertools import combinations

def calculate_tm(primer_sequence):
    primer_sequence = primer_sequence.upper()
    a_t_count = primer_sequence.count('A') + primer_sequence.count('T')
    g_c_count = primer_sequence.count('G') + primer_sequence.count('C')
    return 2 * a_t_count + 4 * g_c_count

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([complement[base] for base in seq[::-1].upper()])

def check_structure(sequence):
    reversed_seq = sequence[::-1]
    hairpin = sequence in reversed_seq or reversed_seq in sequence
    return "Yes" if hairpin else "No"

def advanced_dimer_score(primer1, primer2):
    primer2_rc = reverse_complement(primer2)
    max_score = 0
    for i in range(len(primer1)):
        for j in range(len(primer2_rc)):
            length = 0
            gc_count = 0
            while (i + length < len(primer1)) and (j + length < len(primer2_rc)) and (primer1[i + length] == primer2_rc[j + length]):
                if primer1[i + length] in 'GC':
                    gc_count += 1
                length += 1
            score = length + gc_count
            max_score = max(max_score, score)
    return max_score

def analyze_primers(excel_file_path, output_file_path):
    primers_df = pd.read_excel(excel_file_path)
    primers_df['Fp Tm (°C)'] = primers_df['ForwardPrimer(Fp)'].apply(calculate_tm)
    primers_df['Rp Tm (°C)'] = primers_df['ReversePrimer(Rp)'].apply(calculate_tm)
    primers_df['Fp Hairpin'] = primers_df['ForwardPrimer(Fp)'].apply(check_structure)
    primers_df['Rp Hairpin'] = primers_df['ReversePrimer(Rp)'].apply(check_structure)

    primer_sequences = list(primers_df['ForwardPrimer(Fp)']) + list(primers_df['ReversePrimer(Rp)'])
    dimer_scores = {}
    for primer1, primer2 in combinations(primer_sequences, 2):
        score = advanced_dimer_score(primer1, primer2)
        if score > 0:
            dimer_scores[(primer1, primer2)] = score

    for index, row in primers_df.iterrows():
        fp = row['ForwardPrimer(Fp)']
        rp = row['ReversePrimer(Rp)']
        fp_scores = [score for (p1, p2), score in dimer_scores.items() if p1 == fp or p2 == fp]
        rp_scores = [score for (p1, p2), score in dimer_scores.items() if p1 == rp or p2 == rp]
        primers_df.at[index, 'Max Dimer Score (Fp)'] = max(fp_scores) if fp_scores else 0
        primers_df.at[index, 'Max Dimer Score (Rp)'] = max(rp_scores) if rp_scores else 0

    primers_df.to_excel(output_file_path, index=False)
    return primers_df

if __name__ == "__main__":
    input_file_path = 'pms2_primers(1).xlsx'  # file path
    output_file_path = 'pms2_primers_results(1-2).xlsx'  # desired output file path
    result_df = analyze_primers(input_file_path, output_file_path)
    print(f"Results exported to {output_file_path}")
