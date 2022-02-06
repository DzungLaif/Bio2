from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq

def MSA():
    number = input("Number of protein sequences: ")
    seqrs = []
    for i in range(number):
        temp = input(f'Enter the sequence number {i} :')
        seqr_temp = SeqRecord(Seq(temp), id=str(i))
        seqrs.append(seqr_temp)
    align = MultipleSeqAlignment(*seqrs)
    print(align)
