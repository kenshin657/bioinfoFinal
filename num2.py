import random

from Bio.HMM.MarkovModel import HiddenMarkovModel
from Bio.Seq import Seq
from Bio.Alphabet import *
from  Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO
from pomegranate import *
from pomegranate.distributions.DiscreteDistribution import *
from pomegranate.distributions.DiscreteDistribution import DiscreteDistribution
import pomegranate

#PIP INSTALL BIOPYTHON & PIP INSTALL POMEGRANATE

"""""""""""""""""""""""""""""
PUT FILES IN A SEQUENCE
"""""""""""""""""""""""""""""

records = SeqIO.read("bact.fasta", "fasta")

bacteriaPhage = records.seq

print(bacteriaPhage)

#FOR CG
d0 = DiscreteDistribution({'A':0.2462, 'C':0.2476, 'G':0.2985, 'T':0.2077})
#FOT AT
d1 = DiscreteDistribution({'A':0.2700, 'C':0.2084, 'G':0.1980, 'T':0.3236})

#HMM SETUP
s1 = State(d0, name='CG')
s2 = State(d1, name='AT')


hmm = HiddenMarkovModel('Regions')
hmm.add_states(s1, s2)
hmm.add_transition(hmm.start, s1, 0.5)
hmm.add_transition(hmm.start, s2, 0.5)
hmm.add_transition(s1, s1, 0.9998)
hmm.add_transition(s2,s2, 0.9998)
hmm.add_transition(s1,s2, 0.0002)
hmm.add_transition(s2,s1, 0.0002)
hmm.bake()

viterbi = hmm.predict( bacteriaPhage, algorithm='viterbi')[1:-1]

viterbiPath =''.join(map(str, viterbi))

print(viterbiPath)

