import operator
from math import log
from random import randint
import random

nt_code = {'A':0, 'C':1, 'G':2, 'T':3}

def HammingDistance(Pattern1, Pattern2):
    HDistance = 0
    for i in xrange(0, len(Pattern1)):
        if Pattern1[i] != Pattern2[i]:
            HDistance += 1
    return HDistance


def Neighbors(Pattern, d):
    set = ['A', 'C', 'G', 'T']
    if d == 0:
        return [Pattern]
    if len(Pattern) == 1:
        return ['A', 'C', 'G', 'T']
    Neighborhood = []
    SuffixNeighbors = Neighbors(Pattern[1:], d)
    for text in SuffixNeighbors:
        if HammingDistance(Pattern[1:], text) < d:
            for x in set :
                Neighborhood.append(''.join([x,text]))
        else:
                Neighborhood.append(''.join([Pattern[0],text]))
    return Neighborhood


def MOTIFENUMERATION(Dna, k, d):
        Patterns = []
        for one_sequence in Dna[:]:
            for i in xrange(0, len(one_sequence) - k + 1):
                Pattern = one_sequence[i:i+k]
                neibours = Neighbors(Pattern, d)
                for neibour in neibours:
                    consistant = 0
                    for substring in Dna[:]:
                        check = ['']
                        check = ['yes' for j in xrange(0, len(substring) - k + 1) if HammingDistance(neibour, substring[j:j+k]) <= d]
                        if 'yes' in check:
                            consistant = 1
                        else:
                            consistant = 0
                            break
                    if consistant == 1:
                        Patterns.append(neibour)
        Patterns = list(set(Patterns))
        return Patterns


def MostFreqWord(Column):
    freq = {}
    for element in Column:
        freq[element] = 0
    for element in Column:
        freq[element] += 1
    freq_list = sorted(freq.items(), key=operator.itemgetter(1), reverse=True)
    return freq_list[0][0]


def Score(Motifs):
    Columns = []
    Score = 0
    for i in xrange(0, len(Motifs[0])):
        pattern = []
        for motif in Motifs:
            pattern.append(motif[i])
        Columns.append(pattern)
    for column in Columns:
        consensus_word = MostFreqWord(column)
        for word in column:
            if word != consensus_word:
                Score += 1
    return Score


def Profile(Motifs, Type):
    Columns = []
    for i in xrange(0, len(Motifs[0])):
        pattern = []
        for motif in Motifs:
            word = motif[i].upper()
            pattern.append(word)
        Columns.append(pattern)
    Word_freq = {}
    Word_freq_list = []
    for letter in ['A', 'C', 'G', 'T']:
        Word_freq[letter] = []
    for column in Columns:
        Word_count = {}
        Total_count = 0
        for letter in ['A', 'C', 'G', 'T']:
            Word_count[letter] = 0
        for word in column:
            Word_count[word] += 1
            Total_count += 1
        for letter in ['A', 'C', 'G', 'T']:
            Word_freq[letter].append(float(Word_count[letter]) / float(Total_count))
    for letter in ['A', 'C', 'G', 'T']:        
        Word_freq_list.append(Word_freq[letter])
    if Type == 'list':
        return Word_freq_list
    else:
        return Word_freq


def Entropy(Motifs):
    motif_profile = Profile(Motifs)
    entropy_score = 0.0
    for i in xrange(0, len(motif_profile['A'])):
        column_entropy = 0.0
        for letter in ['A', 'C', 'G', 'T']:
            word_freq = motif_profile[letter][i]
            if word_freq == 0:
                continue
            else:
                column_entropy -= word_freq*log(word_freq, 2)
                entropy_score -= word_freq*log(word_freq, 2)
        print column_entropy
    return entropy_score


def NumberToPattern(index, k):
    nt_anticode = {0:'A',1:'C',2:'G',3:'T'}
    Pattern = []
    for i in xrange(0,k):
        n = (index % 4**(i+1))/4**i
        Pattern.insert(0,nt_anticode[n])
    Pattern_str = ''.join(Pattern)
    return Pattern_str


def DistanceStr(Pattern, Text):
    Distance = 99999
    k = len(Pattern)
    for i in xrange(0, len(Text) - k + 1):
        if HammingDistance(Pattern, Text[i:i+k]) < Distance:
            Distance = HammingDistance(Pattern, Text[i:i+k])
    return Distance


def DistanceList(Pattern, Dna):
    Distance = 0
    for Text in Dna:
        Distance = Distance + DistanceStr(Pattern, Text)
    return Distance


def MEDIANSTRING(Dna, k):
    distance = 99999
    arraySize = 4**k
    for i in xrange(0, arraySize):
        Pattern = NumberToPattern(i, k)
        if distance > DistanceList(Pattern, Dna):
            distance = DistanceList(Pattern, Dna)
            Median = Pattern
    return Median


def ProfileMostProbablekmer(Text, k, Profile):
    candiate_prob = 0.0
    candidate_string = Text[0:k]
    for i in xrange(0, len(Text) - k + 1):
        string = Text[i:i+k]
        propability = 1.0
        for j in xrange(0, k):
            propability = propability * Profile[nt_code[string[j]]][j]
        if propability > candiate_prob:
            candiate_prob = propability
            candidate_string = string
    return candidate_string


def ProfileByLaplacesRule(Motifs, Type):
    Columns = []
    for i in xrange(0, len(Motifs[0])):
        pattern = []
        for motif in Motifs:
            word = motif[i].upper()
            pattern.append(word)
        Columns.append(pattern)
    Word_freq = {}
    Word_freq_list = []
    for letter in ['A', 'C', 'G', 'T']:
        Word_freq[letter] = []
    for column in Columns:
        Word_count = {}
        Total_count = 0
        for letter in ['A', 'C', 'G', 'T']:
            Word_count[letter] = 1   #instead of 0
        for word in column:
            Word_count[word] += 1
            Total_count += 1
        for letter in ['A', 'C', 'G', 'T']:
            Word_freq[letter].append(float(Word_count[letter]) / float(Total_count))
    for letter in ['A', 'C', 'G', 'T']:        
        Word_freq_list.append(Word_freq[letter])
    if Type == 'list':
        return Word_freq_list
    else:
        return Word_freq


def GREEDYMOTIFSEARCH(Dna, k, t):
    BestMotifs = [text[0:k] for text in Dna]
    for i in xrange(0, len(Dna[0]) - k + 1):
        Motif = Dna[0][i:i+k]
        all_motifs = []
        all_motifs.append(Motif)
        for i in xrange(2, t + 1):
            profile = ProfileByLaplacesRule(all_motifs, 'list')
            all_motifs.append(ProfileMostProbablekmer(Dna[i-1], k, profile))
        if Score(all_motifs) < Score(BestMotifs):
            BestMotifs = all_motifs
    return BestMotifs


def MatchMotifs(my_profile, Dna):
    k = len(my_profile[0])
    matched_motifs = []
    for i in xrange(0, len(Dna)):
        best_score = 0.0
        best_motif = ''
        for j in xrange(0, len(Dna[i]) - k + 1):
            motif = Dna[i][j:j+k]
            score = 1.0
            for l in xrange(0, k):
                score = score * my_profile[nt_code[motif[l]]][l]
            if score > best_score:
                best_score = score
                best_motif = motif
        matched_motifs.append(best_motif)
    return matched_motifs


def RANDOMIZEDMOTIFSEARCH(Dna, k, t):
    Motifs = []
    for i in xrange(0, t):
        index  = randint(0,len(Dna[i]) - k - 1)
        Motifs.append(Dna[i][index:index+k])
    BestMotifs = Motifs
    while 1:
        my_profile = ProfileByLaplacesRule(Motifs, 'list')
        Motifs = MatchMotifs(my_profile, Dna)
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs = Motifs
        else:
            return BestMotifs


def roll(massDist):
    randRoll = random.random() # in [0,1)
    sum = 0.0
    result = 1
    for mass in massDist:
        sum += mass
        if randRoll < sum:
            return result
        result += 1


def GIBBSSAMPLER(Dna, k, t, N):
    Motifs = []
    for i in xrange(0, t):
        index  = randint(0,len(Dna[i]) - k - 1)
        Motifs.append(Dna[i][index:index+k])
    BestMotifs = Motifs
    for j in xrange(0, N):
        i = randint(0, t-1)
        if i < t-1:
            new_motifs = Motifs[:i] + Motifs[i+1:]
        else:
            new_motifs = Motifs[:i]
        my_profile = ProfileByLaplacesRule(new_motifs, 'list')
        random_motif = []
        scores = []
        for l in xrange(0, len(Dna[i]) - k + 1):
            motif = Dna[i][l:l+k]
            score = 1.0
            for m in xrange(0, k):
                score = score * my_profile[nt_code[motif[m]]][m]
            scores.append(score)
        scores = [x / sum(scores) for x in scores]
        bd_index  = roll(scores) - 1
        Motifi = Dna[i][bd_index:bd_index+k]
        Motifs = Motifs[:i] + [Motifi] + Motifs[i+1:]
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


input_file = 'input.txt'
with open(input_file, 'r') as f:
    all_motif = []
    for line in f:
        mytext = line.rstrip('\n')
        mytext = mytext.upper()
        mylist = mytext.split(' ')
        #mylist = [float(x) for x in mylist]
        all_motif.append(mytext)
    #print GIBBSSAMPLER(all_motif, 15, 20, 20)
    #print MEDIANSTRING(all_motif, 7)
    best_score = 99999
    best_motifs = []
    for i in xrange(0,200):
        Motifs = GIBBSSAMPLER(all_motif, 7, 14, 1000)
        if Score(Motifs) < best_score:
            best_score = Score(Motifs)
            best_motifs = Motifs
    print best_motifs
    
    #pattern =  MEDIANSTRING(all_motif, 3)
    #print pattern
    #for motif in all_motif:
    #    print DistanceStr('GGGG', motif)
    #print ProfileMostProbablekmer('CGTTGGGTCGAAAAAGCTAAGGAGGATGAAACAGCCCCATATTCTACAACCGCCAACAGTGAGCAATTTCCTACATAGTACTCGGAAACCTGTAATACTATTCTGATTGCCTAGTATGGATGTACCACCAATCCTCTGCCCGTCTCAATAACCATACCGGTCTAGCTATTTGTCTACGGAGGGCGGTTGGTGAGGATTGTATGGAAATTACCTCCATCTAATCTCTGCAGGTTCCTTGACGTAACAAGAATATTATCCACACTACTTTACCCTCTAGTCCAACTGCACCCAGCCCGCGCAGCGTAAGTTCTACCACCCCCTTAATTAGCTCTTGCGCCGCGATAACAGCATTCAGTGAATAGAGTCTTAAGTTAATCGCATACCCTGTAGTGGTGATTCTACTCCCATAACGGCACGACCTCCAAAGTTATTTGCTTTATTTGAGACTATGGACGTTAAACATCTCCCGGAGCGCCGATAAGGTGGATCGTAGCACTACCTCCGCAGATGGTGTGTTAAAGCCCTGCAGACCAACCATTGCGTAGAATATTCTAGCTCTGAAGTTAAACGTACCGTCCGAAAAACTAGGCGCCAGAGGGGGTAACGGATAGTTTAGTAGAATCCATAGAGCGGCGGTTCGATCGACGCCGACAACTGCTGACTAGGCTTACATTGATATTCTACCGTCAAGGGATCCATTAGAGCCTGTGATCGTCCCGACCTAGAGAATAGTTTACAGCATTGAATCCGACTGCCTAGGTAGTGCTTTAACAATGCTGTGTTACACGACTTACGACGCACCAATGCAGGAAACGCAGCGCATTCCATTCGCCCCAGAATGCTCCCTTACAGACACGCGGCTGGTGACTCTGGCCAATGGCAGGAGAACCCAAAACCTTTGATTCGAATTCGGTTTGGTAAAGATCAGCTCGATGTCGGTTAGCATATGGTCACATCTCCGAACTAGCCTTCATGATGGTAACATGTGGCCATACGCTGT',14,all_motif)