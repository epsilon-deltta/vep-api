import random

def get_random_dna():
    return ''.join(random.choices(['a', 'g', 'c', 't'],k=random.randrange(1,15) ) )

def get_percentage(dna):
    dna_percent={'a':0.,'t':0.,'c':0.,'g':0.}
    for s in dna:
        dna_percent[s] += 1./len(dna)*100
                
    return dna_percent