from Bio.Seq import Seq
from Bio import motifs

instances = []
instances.append(Seq("TATAA"))
instances.append(Seq("TATTA"))
instances.append(Seq("TTTAT"))
instances.append(Seq("TATAC"))

m = motifs.create(instances)

print(type(m))
print(m)
print(len(m))
print(m.consensus)
print(m.pwm)
print(m.counts)
print(m.pssm)

m.weblogo("mymotif.png")

pwm = m.counts.normalize(pseudocounts=0.5)
pssm = pwm.log_odds()
print(pwm)
print(pssm)

# exact matches of the instances

test_seq = Seq("TTTTATACACTGCATATAACAACCCAAGCATTATAA")

for pos, seq in m.instances.search(test_seq):
    print(pos, " ", seq)

# using PSSM to score matches
for position, score in pssm.search(test_seq, threshold=4.0):
    print("Position %d: score = %5.3f" % (position, score))

# scores for all positions
print(pssm.calculate(test_seq))
