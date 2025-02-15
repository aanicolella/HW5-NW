# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    alignScores = {}
    # align for each species and add to scores dict
    alignScores['gg'] = nw.align(hs_seq, gg_seq)[0]
    alignScores['mm'] = nw.align(hs_seq, mm_seq)[0]
    alignScores['br'] = nw.align(hs_seq, br_seq)[0]
    alignScores['tt'] = nw.align(hs_seq, tt_seq)[0]

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # sort scores highest to lowest then print in resulting order
    alignScores_sorted = dict(sorted(alignScores.items(), key=lambda x: x[1], reverse=True))
    print('Resulting alignment scores:')
    for species, score in alignScores_sorted.items():
        print(f'{species} score: {score}')
    

if __name__ == "__main__":
    main()
