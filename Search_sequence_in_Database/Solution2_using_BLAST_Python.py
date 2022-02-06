from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

def blast_nr_prot(filename):
    record = SeqIO.read(open(filename), format="fasta")
    result = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))
    save_file = open("analyse_result.xml", "w")
    save_file.write(result.read())
    save_file.close()
    result.close()

# use blast record to alanyse the result from the function above
# input file type: xml
def blast_records(xml_filename):
    blastRecords = NCBIXML.read(open(xml_filename))
    E_threshold = 0.001
    for blastRecord in blastRecords:
        for alignment in blastRecord.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_threshold:
                    print("∗∗∗∗Alignment∗∗∗∗")
                    print("sequence:", alignment.title)
                    print("length:", alignment.length)
                    print("e value:", hsp.expect)
                    print(hsp.query[0:75] + "...")
                    print(hsp.match[0:75] + "...")
                    print(hsp.sbjct[0:75] + "...")

    # Check the first alignment
    first_alignment = blastRecords.alignments[0]
    print("FIRST ALIGNMENT:")
    print("Accession: " + first_alignment.accession)
    print("Hit id: " + first_alignment.hit_id)
    print("Definition: " + first_alignment.hit_def)
    print("Alignment length: ", first_alignment.length)
    print("Number of HSPs: ", len(first_alignment.hsps))

    # Check single HSP
    hsp = first_alignment.hsps[0]
    print("E−value: ", hsp.expect)
    print("Score: ", hsp.score)
    print("Length: ", hsp.align_length)
    print("Identities: ", hsp.identities)
    print("Alignment of the HSP:")
    print(hsp.query)
    print(hsp.match)
    print(hsp.sbjct)

    # print top 10 alignment
    print("Top 10 alignments:")
    for i in range(10):
        alignment = blastRecords.alignments[i]
    print("Accession: " + alignment.accession)
    print("Definition: " + alignment.hit_def)
    for hsp in alignment.hsps:
        print("E−value: ", hsp.expect)
    print()

    # Find which organisms in top 20 alignments
    import re
    specs = []
    for i in range(20):
        alignment = blastRecords.alignments[i]
        definition = alignment.hit_def
        x = re.search("\[(.∗?)\]", definition).group(1)
        specs.append(x)

    print("Organisms: ")
    for s in specs: print(s)