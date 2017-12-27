from itertools import product

bases = 'ACGT'

mutations=[]
for trimer in product(bases, bases, bases):
    for base in bases:
        if trimer[1] != base:
            mutations.append((''.join(trimer), base))
