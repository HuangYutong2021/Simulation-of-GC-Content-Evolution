import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import statsmodels.api as sm
from collections import deque
import random
import scipy.stats

Corresponding={'A':'T','T':'A','C':'G','G':'C'}

# Gram-negative
package_candidates=[
            'Helicobacteraceae_bacterium(SAMEA5278417)',
            'Candidatus_Helicobacter_avistercoris(SAMN15816903)',
# ... up to your named files
            ]
filename_candidates=[
                     '2212472.SAMEA5278417.genes',
                     '2838619.SAMN15816903.genes',
# ... up to your named files
            ]


def read_data(index):
    package=package_candidates[index]
    file_name=filename_candidates[index]
    number=file_name.split(".")[0]
    species_name=package.split("(")[0].replace("_"," ")
    gene_file=f"../{package}/{file_name}.fasta"
    disorder_file=f'../{package}/{number}.result'
    print(package)
    Total_seq = []
    with open(gene_file) as fa:
        seq = ''
        for line in fa:
            if not line.startswith('>'):
                seq += line.replace('\n', '')
            elif len(seq) != 0:
                if len(seq) % 3 == 0:
                    Total_seq.append(seq)
                seq = ''
    Total_seq.append(seq)
    return Total_seq,disorder_file, species_name


def Count_disorder_per_residue(disorder_file):
    Total_dis_n = 0
    Total_count = 0
    with open(disorder_file) as f:
        for line in f:
            if line[0].isdigit():
                value=float(list(line.split("\t"))[2])
                if value > 0.5:
                    Total_dis_n+=1
                Total_count+=1
    return Total_dis_n,Total_count

def Count_disorder_per_window(disorder_file):
    Total_dis_n = 0
    Total_count = 0
    Len_limit = 15
    chunk_count = 0
    chunk_dis = 0
    with open(disorder_file) as f:
        for line in f:
            if line[0].isdigit():
                if int(list(line.split("\t"))[0]) == 1:
                    Total_count += chunk_count
                    Total_dis_n += chunk_dis
                    chunk_count = 0
                    chunk_dis = 0
                    window = deque(maxlen=Len_limit)
                value = float(list(line.split("\t"))[2])
                window.append(value)
                if len(window) == Len_limit:
                    chunk_count += 1
                    if np.mean(window) > 0.5:
                        chunk_dis += 1
    return Total_dis_n, Total_count



def Count_CG_percentage(Total_seq): # Total_seq=['xxxx','xxxx']
    GC_n=0
    Total_n=0
    for seq in Total_seq:
        for i in seq:
            if i=='C' or i=='G':
                GC_n+=1
            Total_n+=1
    return GC_n/Total_n


### FORMALLY BEGIN ###

Codens={
  'TTT':'F','TTC':"F",'TTA':'L','TTG':"L",
  'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
  'TAT':'Y','TAC':'Y','TAA':'STOP','TAG':'STOP',
  'TGT':'C','TGC':'C','TGA':'STOP','TGG':'W',
  'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
  'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
  'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
  'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
  'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
  'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
  'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
  'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
  'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
  'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
  'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
  'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
  'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
  'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
  'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
  'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
Dict={} # records aa:potential codens
for key, value in Codens.items():
    if value not in Dict.keys():
        Dict[value]=[key]
    else:
        Dict[value].append(key)

'''
Change coden to see the consequences
'''
def create_mutation(prob,total_seq): # prob: The probability of create mutation for each amino acid
    new_total=[]
    Count=0
    for seq in total_seq:
        Count+=1
        new_seq=''
        for i in range(0, len(seq), 3):
          try:
              if Codens[seq[i:i + 3]] != 'STOP':
                  aa = Codens[seq[i:i + 3]]
                  potential = list(set(Dict[aa]) - set([seq[i:i + 3]]))  # Except itself
                  if len(potential) != 0:
                      if np.random.choice([0, 1], 1, p=[1 - prob, prob])[0] == 1:
                          P = np.array([1 / len(potential) for _ in range(len(potential))])  # 剩下的codens 概率均分
                          new_seq += np.random.choice(potential, p=P.ravel())
                      else:
                          new_seq += seq[i:i + 3]
                  else:
                      new_seq += seq[i:i + 3]
          except KeyError:
              pass
        new_total.append(new_seq)
    return new_total


def create_mutation_mode3(total_seq): # prob: The probability of create mutation for each amino acid
    new_total=[]
    Count=0
    for seq in total_seq:
        Count+=1
        new_seq=''
        for i in range(0, len(seq), 3):
          try:
              if Codens[seq[i:i + 3]] != 'STOP':
                  aa = Codens[seq[i:i + 3]]
                  potential=Dict[aa] # All neutral mutation (include itself)
                  P = np.array([1 / len(potential) for _ in range(len(potential))])  # Equipartition of probability
                  new_seq += np.random.choice(potential, p=P.ravel())

          except KeyError:
              pass
        new_total.append(new_seq)
    return new_total



def create_mutation_extreme(ishighest,total_seq): # prob: The probability of create mutation for each amino acid
    new_total=[]
    Count=0
    for seq in total_seq:
        Count+=1
        new_seq=''
        for i in range(0, len(seq), 3):
          try:
            if Codens[seq[i:i + 3]] != 'STOP':
                aa= Codens[seq[i:i + 3]]
                ordered_list=sorted(Dict[aa], key=lambda s: s.count("C")+s.count("G"))
                if ishighest:
                    new_seq+=ordered_list[-1]
                else:
                    new_seq += ordered_list[0]
          except KeyError:
              pass
        new_total.append(new_seq)
    return new_total



def create_mutation_extreme_given(case_choice,total_seq): # prob: The probability of create mutation for each amino acid
    new_total=[]
    Count=0
    for seq in total_seq:
        Count+=1
        new_seq=''
        for i in range(0, len(seq), 3):
          try:
            if Codens[seq[i:i + 3]] != 'STOP':
                aa= Codens[seq[i:i + 3]]
                ordered_list=sorted(Dict[aa], key=lambda s: s.count("C")+s.count("G"))
                if case_choice:
                    new_seq+=ordered_list[-1]
                else:
                    new_seq += ordered_list[0]
          except KeyError:
              pass
        new_total.append(new_seq)
    return new_total


def main():
    GC_all=[]
    disorder_all=[]
    species_names=[]
    for i in range(len(package_candidates)):
        Total_seq, disorder_file,species_name = read_data(i)
        species_names.append(species_name)
        GC_content=Count_CG_percentage(Total_seq)
        GC_all.append(GC_content)
        print('GC content_percentage',GC_content)
        Total_dis_n,Total_count=Count_disorder_per_residue(disorder_file)
        print('disorder_percentage',Total_dis_n/Total_count)
        disorder_all.append(Total_dis_n/Total_count)


    X = np.array(disorder_all).reshape(-1, 1) # converted to 1 column which is a requirement
    print(X)
    y = np.array(GC_all)

    # Apply linear regression
    X1 = sm.add_constant(X)
    model = sm.OLS(y, X1).fit()
    print('model summary')
    print(model.summary())

    p_values = model.pvalues
    print("P-values of the model parameters:")
    print("p-value",p_values)
    print("model parameters:")
    print(model.params[0],model.params[1])

    # Find outliers
    outliers = np.where(np.abs(standardized_residuals) > 2)
    print("Standard residuals",outliers)


    # Calculate Cook's D value
    influence = model.get_influence()
    cd, _ = influence.cooks_distance
    df = pd.DataFrame(list(zip(species_names,cd)), columns = ['species',"Cook's D"])
    print(df)

    n=len(X)
    # k = model.df_model
    # threshold = 4 / (n - k - 1)
    threshold = 4 / n
    print("threshold",threshold)
    outliers = np.where(cd >= threshold)
    # Show the outliers
    print("outliers",outliers)

    # LinearRegression
    model = LinearRegression()
    model.fit(X, y)
    r_squared = model.score(X, y)
    print("original R squared",r_squared)
    print("pearson",scipy.stats.pearsonr(disorder_all, y)[0])
    return disorder_all




def main_mutated(case_choice):
    GC_all = []
    disorder_all = []
    species_names = []
    for i in range(len(package_candidates)):
        Total_seq, disorder_file, species_name = read_data(i)
        species_names.append(species_name)
        if case_choice=="high":
            new_seq = create_mutation_extreme(True, Total_seq)
        elif case_choice=="low":
            new_seq = create_mutation_extreme(False, Total_seq)
        elif case_choice=="Mode3":
            new_seq = create_mutation_mode3(Total_seq)
        elif type(case_choice)==list:
            if len(case_choice)==len(package_candidates):
                new_seq = create_mutation_extreme_given(case_choice[i],Total_seq)
            else:
                raise ValueError("Please recheck your case_choice, maybe the length is wrong")
        else:
            new_seq = create_mutation(case_choice, Total_seq)


        GC_content = Count_CG_percentage(new_seq)
        print(f'{species_name} GC content', GC_content)
        GC_all.append(GC_content)
        Total_dis_n, Total_count = Count_disorder_per_residue(disorder_file)
        print('disorder_percentage', Total_dis_n / Total_count)
        disorder_all.append(Total_dis_n / Total_count)

    X = np.array(disorder_all).reshape(-1, 1)  # converted to 1 column which is a requirement
    print(X)
    y = np.array(GC_all)

    # Linear Regression
    X1 = sm.add_constant(X)
    model = sm.OLS(y, X1).fit()
    print('model summary')
    print(model.summary())

### USER FACE ###

# Analyze the original data, and investigate its relationship.
disorder_all=main()
'''
choice: 'high': Extreme highest GC,
        'low': Extreme lowest GC,
        'Mode3':Treat all codens equally,
        0-1 any float number (normal randomization mutation probability),
        or given choice like choice=[0,1,1,0,0,1,1,0,0,0,1,0,0,1,0]
'''
# Determine your choice
choice="Mode3"
# Make mutations!!! and you can investigate its relationship.
main_mutated(choice)