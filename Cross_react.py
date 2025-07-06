import random
import pandas as pd
import numpy as np
import subprocess
from pybatman.functions import train, peptide2index
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
from abc import ABC, abstractmethod
from mhctools import MHCflurry

class Population:
    """A class to represent a population.
    Attributes:
    path_to_IEPAPI: the absolute path to IEPAPI folder
    freq: the absolute path to tsv table with genotype frequencies across population. Tsv file consist of 2 columns without headers. First is a column with absolute counts of each genotype, second is genotypes themselfs
    path_to_hla_db: absolute path to table with hla pseudosequences. The format of the table is the same as in IEPAPI table (data/pseudoSequence(ELIM).csv file) """
    
    def __init__(self,freq, path_to_IEPAPI=None, path_to_hla_db=None):
        if path_to_IEPAPI is None:
            path_to_IEPAPI = os.getcwd() + '/IEPAPI/'
        if path_to_hla_db is None:
            self.IEPAPI_db = pd.read_csv(path_to_IEPAPI + 'data/pseudoSequence(ELIM).csv' )
        else:
            self.IEPAPI_db = pd.read_csv(path_to_hla_db)
        if "tsv" in freq:
            self.genotype_freq = self.pars_tsv(freq)
        else:
            self.genotype_freq = freq
        new_genotype_df = self.genotype_freq
        c = 0
        for i in range(len(self.genotype_freq)):
            genotype = self.genotype_freq['Genotype'][i]
            if any(hla not in self.IEPAPI_db['HLA'].tolist() for hla in [('HLA-'+hla[0:7]).replace('*','').replace(':','') for hla in genotype]):
                c += self.genotype_freq['Count'][i]
                new_genotype_df = new_genotype_df.drop(i)
                print(f"Warning! The genotype {genotype} was removed because it contains unknown alleles")
        new_genotype_df['Freq'] = np.array(new_genotype_df['Freq'])*sum(self.genotype_freq['Count'])/(sum(self.genotype_freq['Count'])-c)
        self.genotype_freq = new_genotype_df.reset_index(drop=True)
    
    
    def pars_tsv(self,path_to_tsv):
        df = pd.read_csv(path_to_tsv,sep='\t',names=['ID','HLA'])
        ids = np.unique(df['ID'])
        freq_d = {'Genotype':[],'Count':[],'Freq':None}
        for i in ids:
            genotype = df.loc[df['ID'] == i]['HLA'].tolist()
            genotype.sort()
            if genotype in freq_d['Genotype']:
                freq_d['Count'][freq_d['Genotype'].index(genotype)] += 1
            else:
                freq_d['Genotype'].append(genotype)
                freq_d['Count'].append(1)
        freq_d['Freq'] = np.array(freq_d['Count'])/len(ids)
        return pd.DataFrame(freq_d)
    
    
    def get_individual_sample(self,immun_predictor,N=1,mean_number=10):
        sample_list = []
        for n in range(N):
            genotype = np.random.choice(a=np.array(self.genotype_freq['Genotype']), p = self.genotype_freq['Freq'])
            sample_list.append(Individual(genotype,immun_predictor,mean_number))
        return sample_list



class Protein:
    """A class to represent a protein
    Attributes:
    protein_seq: amino acid sequence of protein"""
    
    def __init__(self,seq):
        self.protein_seq = seq
        self.peptide_pool = None



class Protein_processing:

    def get_peptide_weight(self):
        return None
    
    def processing(self, protein, k=9):
        peptide_weight = self.get_peptide_weight()
        peptide_list = []
        seq = protein.protein_seq
        for i in range(len(seq)-k):
            peptide_list.append(seq[i:i+k])
        if peptide_weight is not None:
            peptide_list = np.random.choice(a = peptide_list,p = peptide_weight)
        protein.peptide_pool = sorted(peptide_list)
        return protein


class peptide_MHC_complex:
    """A class to represent pMHC complex
    Attributes:
    peptide: peptide sequence
    HLA: HLA name
    pres: presentation score for peptide-HLA pair
    immun: immunogenecity score for peptide-HLA pair
    HLA_i: HLA index"""
    
    def __init__(self,peptide,HLA,pres,immun,HLA_i):
        self.peptide = peptide
        self.HLA = HLA
        self.pres = pres
        self.immun = immun
        self.HLA_i = HLA_i
    
    def cross (self,cross_react_predictor,mutant_pMHC):
        homology = cross_react_predictor.predict([self.peptide],[mutant_pMHC.peptide])
        c=0
        for i in range(len(self.peptide)):
            if self.peptide[i] != mutant_pMHC.peptide[i]:
                c += 1
        if c > 2:
            homology = 0
        homology = homology*mutant_pMHC.pres/self.pres
        if homology > 1:
            homology = 1
        return homology


class Individual:
    """A class to represent a individual
    Attributes:
    genotype: list with HLA alleles
    immun_predictor: object of the class IEPAPI_immunogenicity_predictor
    mean_number: mean number of pMHC complexes that will be presented
    protein_wild: object of the class Protein that represented a wild protein
    """
    
    def __init__(self,genotype,immun_predictor,mean_number):
        self.genotype = sorted(genotype)
        self.protein_wild = None
        
    def get_pMHC_pool(self,protein, immun_predictor):
        pep_list = protein.peptide_pool
        pMHC_list = []
        pMHC_d = {'HLA':[],'peptide':[],'pres':[],'immun':[],'HLA_i':[]}
        for j,hla in enumerate(self.genotype):
            pMHC_d['HLA_i'] += [j]*len(pep_list)
            pMHC_d['HLA'] += [hla]*len(pep_list)
            pMHC_d['peptide'] += pep_list
        res = immun_predictor.predict(pep_list,self.genotype)
        pMHC_d['pres'] = res[1]
        pMHC_d['immun'] = res[0]
        for i in range(len(pMHC_d['peptide'])):
            pMHC_list.append(peptide_MHC_complex(pMHC_d['peptide'][i],pMHC_d['HLA'][i],pMHC_d['pres'][i],pMHC_d['immun'][i],pMHC_d['HLA_i'][i]))
        return pMHC_list
                
    def set_protein_wild(self, protein_wild,immun_predictor,mean_number=10):
        self.protein_wild = protein_wild
        pool_result = self.get_pMHC_pool(self.protein_wild,immun_predictor)
        self.pMHC_pool = self.pMHC_selection(pool_result,mean_number)
        
    def pMHC_selection(self,pMHC,mean_number):
        P=[]
        for pMHC_i in pMHC:
            P.append(pMHC_i.pres*pMHC_i.immun)
        P = np.array(P)
        if mean_number is not None:
            P = P-((sum(P)-mean_number)/len(P))
            P[P<0] = 0
            P[P>1] = 1
        pMHC_selected = []
        for i in range(len(P)):
            if np.random.choice(a=[True,False], p=[P[i],1-P[i]]):
                pMHC_selected.append(pMHC[i])
        return pMHC_selected
        

    def cross_react(self,cross_react_predictor,immun_predictor,mutant_protein):
        result = []
        wild_pMHC = self.pMHC_pool
        if len(wild_pMHC) == 0:
            return 0
        mutant_pMHC = self.get_pMHC_pool(mutant_protein,immun_predictor)
        for w_pMHC in wild_pMHC:
            for m_pMHC in mutant_pMHC:
                if w_pMHC.HLA_i == m_pMHC.HLA_i:
                    result = np.append(result, w_pMHC.cross(cross_react_predictor,m_pMHC))
        return sum(result)/len(wild_pMHC)


class Basic_cross_react_predictor:

    @abstractmethod
    def predict(self,wild_peptides,mutant_peptides):
        pass

    @abstractmethod
    def norm_distance(self,peptide_distance):
        return peptide_distance

    
class BATMAN_predictor(Basic_cross_react_predictor):
    """A class to BATMAN tool
    Attributes:
    path_to_BATMAN: the absolute path to BATMAN folder (pybatman package)"""
    inferred_AA_matrix, weight_profile, path_to_BATMAN = None,None,None
    max_value = 4.382771069923114
    
    def __init__(self, path_to_BATMAN=None, path_to_train_data=None):
        
        if path_to_BATMAN is None:
            result = subprocess.run(
            ['conda', 'list', '-n', 'BATMAN'],
            capture_output=True,
            text=True,
            check=True)
            for line in result.stdout.splitlines():
                if 'pybatman' in line:
                    parts = line.split()
                    package_path = os.path.join(result.stdout.splitlines()[0].split()[5][:-1], 'lib/python3.11/site-packages/',parts[0])
                    self.path_to_BATMAN = package_path
        else:
            self.path_to_BATMAN = path_to_BATMAN

        if path_to_train_data is None:
            path_to_train_data = os.getcwd() + '/test_input_BATMAN'
            
        full_peptide_data = pd.read_csv(path_to_train_data)
        full_peptide_data.tcr = 'TCR1'
        full_peptide_data.to_csv('BATMAN_train.csv')
#AA_matrix_name = 'blosum100' 
        AA_matrix_prior = pd.read_csv(self.path_to_BATMAN+ '/data/AA_matrices/blosum100.csv',index_col=0)
        inferred_weights, inferred_AA_matrix = train('BATMAN_train.csv','full', AA_matrix_prior, steps=80000, seed=10)
        weight_profile = inferred_weights.to_numpy()
        self.weight_profile = weight_profile
        self.inferred_AA_matrix = inferred_AA_matrix
    
    def predict(self,wild_peptides,mutant_peptides):
        peptide_distance = peptide2index(wild_peptides,
                                mutant_peptides,
                                 self.inferred_AA_matrix,
                                 self.weight_profile)
        peptide_distance = np.array(peptide_distance)
        return self.norm_distance(peptide_distance)


    def norm_distance(self,peptide_distance):
        peptide_distance[peptide_distance > self.max_value] = self.max_value
        peptide_distance[peptide_distance < 0] = 0
        peptide_distance = peptide_distance/self.max_value
        peptide_distance = np.ones(len(peptide_distance)) - peptide_distance
        return peptide_distance



class Basic_immunogenecity_predictor:
    log = {}
    
    @abstractmethod
    def predict(self,peptides,hla_vector):
        pass

    def update_log(self,new_data):
        self.log.update(new_data)

    def get_log(self):
        return self.log

        
class IEPAPI_immunogenicity_predictor(Basic_immunogenecity_predictor):
    """A class to IEPAPI tool
    Attributes:
    path_to_IEPAPI: the absolute path to IEPAPI folder
    python_path: the absolute path to python in IEPAPI environment
    hla_pseudoseq_data: the absolute path to HLA pseudosequence table"""
    hla_pseudoseq_data = None
    path_to_IEPAPI  = None
    python_path = None
    log = {}
    def __init__(self, path_to_IEPAPI=None, path_to_hla_db=None, python_path=None):
        if path_to_IEPAPI is None:
            self.path_to_IEPAPI = os.getcwd() + '/IEPAPI/'
        else:
            self.path_to_IEPAPI = path_to_IEPAPI
        if python_path is None:
            result = subprocess.run(
    ['conda', 'list', '-n', 'IEPAPI'],
            capture_output=True,
            text=True,
            check=True)
            self.python_path = result.stdout.splitlines()[0].split()[5][:-1]+'/bin/python3.7'
        else:
            self.python_path = python_path
        if path_to_hla_db is None:
            self.IEPAPI_db = self.path_to_IEPAPI+'data/pseudoSequence(ELIM).csv'
        else:
            self.IEPAPI_db = path_to_hla_db
        self.hla_pseudoseq_data = pd.read_csv(self.IEPAPI_db)

            
    def predict(self,peptides,hla_vector):
        immun = np.array([-1.0]*(len(peptides)*len(hla_vector)))
        pres = np.array([-1.0]*(len(peptides)*len(hla_vector)))
        with open(self.path_to_IEPAPI+f'IEPAPI_input{os.getpid()}', 'w') as file:
            file.write('peptide,HLA,seq')
            file.write('\n')
            for j,hla in enumerate(hla_vector):
                if 'HLA-' not in hla:
                    hla = 'HLA-'+hla
                if len(hla) > 11:
                    hla = hla[0:11]
                hla_to_df = hla.replace('*','').replace(':','')
                hla_pseudoseq = self.hla_pseudoseq_data.loc[self.hla_pseudoseq_data['HLA'] == hla_to_df,'pseudoSeq'].iloc[0]
                for i in range(len(peptides)):
                    if peptides[i]+hla in self.log:
                        log_pep = self.log[peptides[i]+hla]
                        pres[i+len(peptides)*j] = log_pep[0]
                        immun[i+len(peptides)*j] = log_pep[1]
                    else:
                        file.write(f'{peptides[i]},{hla},{hla_pseudoseq}')
                        file.write('\n')
        if all(immun != -1):
            os.remove(self.path_to_IEPAPI+f'IEPAPI_input{os.getpid()}')
            return immun, pres
        subprocess.run([self.python_path, self.path_to_IEPAPI+'IEPAPI_predict.py', '--input', f'IEPAPI_input{os.getpid()}', '--output',f'IEPAPI_output{os.getpid()}'],  cwd=self.path_to_IEPAPI)
        out = pd.read_csv(self.path_to_IEPAPI+f'IEPAPI_output{os.getpid()}')
        for j,i in enumerate(np.where(immun==-1.0)[0]):
            pe, h, pr, im = out.iloc[j,0], out.iloc[j,1], out.iloc[j,2], out.iloc[j,3]
            immun[i] = im
            pres[i] = pr
            self.log[pe+h] = (pr,im)
        os.remove(self.path_to_IEPAPI+f'IEPAPI_output{os.getpid()}')
        os.remove(self.path_to_IEPAPI+f'IEPAPI_input{os.getpid()}')
        return immun, pres

    

class Baseline_IEPAPI_predictor(Basic_immunogenecity_predictor):
    """A class to IEPAPI presentation predictor tool
    Attributes:
    path_to_IEPAPI: the absolute path to IEPAPI folder
    python_path: the absolute path to python in IEPAPI environment
    hla_pseudoseq_data: the absolute path to HLA pseudosequence table"""
    hla_pseudoseq_data = None
    path_to_IEPAPI  = None
    python_path = None
    log = {}
    def __init__(self, path_to_IEPAPI=None, path_to_hla_db=None, python_path=None):
        if path_to_IEPAPI is None:
            self.path_to_IEPAPI = os.getcwd() + '/IEPAPI/'
        else:
            self.path_to_IEPAPI = path_to_IEPAPI
        if python_path is None:
            result = subprocess.run(
    ['conda', 'list', '-n', 'IEPAPI'],
            capture_output=True,
            text=True,
            check=True)
            self.python_path = result.stdout.splitlines()[0].split()[5][:-1]+'/bin/python3.7'
        else:
            self.python_path = python_path
        if path_to_hla_db is None:
            self.IEPAPI_db = self.path_to_IEPAPI+'data/pseudoSequence(ELIM).csv'
        else:
            self.IEPAPI_db = path_to_hla_db
        self.hla_pseudoseq_data = pd.read_csv(self.IEPAPI_db)

            
    def predict(self,peptides,hla_vector):
        immun = np.ones(len(hla_vector)*len(peptides))
        pres = np.array([-1.0]*(len(peptides)*len(hla_vector)))
        with open(self.path_to_IEPAPI+f'IEPAPI_input{os.getpid()}', 'w') as file:
            file.write('peptide,HLA,seq')
            file.write('\n')
            for j,hla in enumerate(hla_vector):
                if 'HLA-' not in hla:
                    hla = 'HLA-'+hla
                if len(hla) > 11:
                    hla = hla[0:11]
                hla_to_df = hla.replace('*','').replace(':','')
                hla_pseudoseq = self.hla_pseudoseq_data.loc[self.hla_pseudoseq_data['HLA'] == hla_to_df,'pseudoSeq'].iloc[0]
                for i in range(len(peptides)):
                    if peptides[i]+hla in self.log:
                        log_pep = self.log[peptides[i]+hla]
                        pres[i+len(peptides)*j] = log_pep
                    else:
                        file.write(f'{peptides[i]},{hla},{hla_pseudoseq}')
                        file.write('\n')
        if all(pres != -1):
            os.remove(self.path_to_IEPAPI+f'IEPAPI_input{os.getpid()}')
            return immun, pres
        subprocess.run([self.python_path, self.path_to_IEPAPI+'IEPAPI_predict.py', '--input', f'IEPAPI_input{os.getpid()}', '--output',f'IEPAPI_output{os.getpid()}'],  cwd=self.path_to_IEPAPI)
        out = pd.read_csv(self.path_to_IEPAPI+f'IEPAPI_output{os.getpid()}')
        for j,i in enumerate(np.where(pres==-1.0)[0]):
            pe, h, pr = out.iloc[j,0], out.iloc[j,1], out.iloc[j,2]
            pres[i] = pr
            self.log[pe+h] = pr
        os.remove(self.path_to_IEPAPI+f'IEPAPI_output{os.getpid()}')
        os.remove(self.path_to_IEPAPI+f'IEPAPI_input{os.getpid()}')
        return immun, pres

        
class Baseline_MHCflurry_predictor(Basic_immunogenecity_predictor):
    log = {}
    def predict(self,peptides,hla_vector):
        pres =  np.array([-1.0]*(len(peptides)*len(hla_vector)))
        immun = np.ones(len(hla_vector)*len(peptides))
        pep_to_predict = []
        hla_to_predict = []
        for j,hla in enumerate(hla_vector):
            if 'HLA-' not in hla:
                hla = 'HLA-'+hla
            if len(hla) > 11:
                hla = hla[0:11]
            for i in range(len(peptides)):
                if peptides[i]+hla in self.log:
                    log_pep = self.log[peptides[i]+hla]
                    pres[i+len(peptides)*j] = log_pep
                else:
                    pep_to_predict.append(peptides[i])
                    hla_to_predict.append(hla)
        if all(pres != -1):
            return immun, pres
        out = MHCflurry(alleles=hla_to_predict).predict_peptides(pep_to_predict).to_dataframe()
        idx = np.lexsort((out['peptide'].values, out['allele'].values))
       # out = out.sort_values(by = ['allele','peptide'])
        h_v, pe_v, pr_v = out['allele'].iloc[idx].values, out['peptide'].iloc[idx].values, out['score'].iloc[idx].values
        for j,i in enumerate(np.where(pres==-1.0)[0]):
            h, pe, pr = h_v[j], pe_v[j], pr_v[j]
            pres[i] = pr
            self.log[pe+h] = pr
            
        return immun, pres

class Baseline_homology_predictor(Basic_immunogenecity_predictor):
    log = {}

    def predict(self, peptides, hla_vector):
        immun = np.ones(len(hla_vector*len(peptides)))
        pres = np.ones(len(hla_vector*len(peptides)))
        return immun, pres


def compare_dist(a,b):
    X_samples = np.random.choice(a,1000)
    Y_samples = np.random.choice(b, 1000)
    return np.mean(X_samples > Y_samples)



def ind_dist(individual, cross_react_predictor,immun_predictor,protein1,protein2,k, N_iterations):
    result = []
    for i in range(N_iterations):
        individual.set_protein_wild(protein1, immun_predictor,k)
        result.append(individual.cross_react(cross_react_predictor, immun_predictor,protein2))
    return (result, immun_predictor.get_log())
    
def Cross_react_predict(protein1, protein2, genotypes, immun_predictor,cross_react_predictor,N_iterations=50,k=10,N_proc=2,N_individuals=None):
    result_log = {}
    if type(genotypes).__name__ == 'Population':
        individuals = genotypes.get_individual_sample(immun_predictor,N=N_individuals,mean_number=k)
    else:
        individuals = []
        for genotype in genotypes:
            individuals.append((Individual(genotype,immun_predictor,k)))
    with Pool(processes=N_proc) as pool:
        args = [(ind,cross_react_predictor, immun_predictor, protein1, protein2,k,N_iterations) for ind in individuals]
        result = pool.starmap(ind_dist,args)
    for i in result:
        result_log.update(i[1])
    immun_predictor.update_log(result_log)
    
    return [r[0] for r in result]




def Cross_react_plot(genotypes,immun_predictor,cross_react_predictor, protein1, protein2, protein3=None,protein4=None, protein5=None,N_iterations=50,k=10,N_proc=2,N_individuals=None, savefile=None):
    if N_individuals is None:
        x_len = len(genotypes)
    else:
        x_len = N_individuals
    if protein5 is not None:
        x = np.array([i for i in range(1000, 1000*x_len+1000,1000)])
        c1 = "red"
        c2 = 'blue'
        c3 = 'green'
        c4 = 'purple'
        data1 = Cross_react_predict(protein1, protein2, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        data2 = Cross_react_predict(protein1, protein3, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        data3 = Cross_react_predict(protein1, protein4, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        data4 = Cross_react_predict(protein1, protein5, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        plt.figure()
        plt.boxplot(data1,positions=x-100,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='pink', color=c1),
            capprops=dict(color=c1),
            whiskerprops=dict(color=c1),
            flierprops=dict(color=c1, markeredgecolor=c1),
            medianprops=dict(color=c1))
        plt.boxplot(data2,positions=x+100,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='cyan', color=c2),
            capprops=dict(color=c2),
            whiskerprops=dict(color=c2),
            flierprops=dict(color=c2, markeredgecolor=c2),
            medianprops=dict(color=c2))
        plt.boxplot(data3,positions=x+300,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='springgreen', color=c3),
            capprops=dict(color=c3),
            whiskerprops=dict(color=c3),
            flierprops=dict(color=c3, markeredgecolor=c3),
            medianprops=dict(color=c3))
        plt.boxplot(data4,positions=x-300,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='magenta', color=c4),
            capprops=dict(color=c4),
            whiskerprops=dict(color=c4),
            flierprops=dict(color=c4, markeredgecolor=c4),
            medianprops=dict(color=c4))


#for j,i in enumerate(x):
#    plt.hlines(y = result_baseline_mp_ac[j], xmin=i+25,xmax=i+175, color = 'black')
#    plt.hlines(y = result_baseline_mp_ab[j], xmin=i-175,xmax=i-25, color = 'black')
        plt.xlim(100,1000*x_len+1000)
        plt.title('Cross-reactivity boxplots')
        plt.xticks(x,[str(i) for i in range(x_len)])
        plt.xlabel('Genotype')
        plt.ylabel('Cross-reactivity')
        blue_patch = mpatches.Patch(color='blue', label='C protein')
        red_patch = mpatches.Patch(color='red', label='B protein')
        green_patch = mpatches.Patch(color = 'green',label='D protein')
        purple_patch = mpatches.Patch(color = 'purple',label='E protein')
        plt.legend(handles=[blue_patch, red_patch, green_patch, purple_patch])
#plt.hlines(y = 0.9867, xmin=100,xmax=6000, color = 'red', linestyles='--')
#plt.hlines(y = 0.9867, xmin=100,xmax=6000, color = 'blue', linestyles='--')
#plt.hlines(y = result_base_ab_v, xmin=250,xmax=6000, color = 'red', linestyles='--')
#plt.hlines(y = result_base_ac_v, xmin=250,xmax=6000, color = 'blue', linestyles='--')
#plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,],['5','10','15','20','25','30','40','50','60','70','80','90','100'])
    #plt.savefig('valid_pair.png')
        if savefile is not None:
            plt.savefig(savefile)
            return
        plt.show()
    elif protein4 is not None:
        x = np.array([i for i in range(800, 800*x_len+800,800)])
        c1 = "red"
        c2 = 'blue'
        c3 = 'green'
        c4 = 'purple'
        data1 = Cross_react_predict(protein1, protein2, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        data2 = Cross_react_predict(protein1, protein3, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        data3 = Cross_react_predict(protein1, protein4, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        plt.figure()
        plt.boxplot(data1,positions=x-200,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='pink', color=c1),
            capprops=dict(color=c1),
            whiskerprops=dict(color=c1),
            flierprops=dict(color=c1, markeredgecolor=c1),
            medianprops=dict(color=c1))
        plt.boxplot(data2,positions=x,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='cyan', color=c2),
            capprops=dict(color=c2),
            whiskerprops=dict(color=c2),
            flierprops=dict(color=c2, markeredgecolor=c2),
            medianprops=dict(color=c2))
        plt.boxplot(data3,positions=x+200,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='springgreen', color=c3),
            capprops=dict(color=c3),
            whiskerprops=dict(color=c3),
            flierprops=dict(color=c3, markeredgecolor=c3),
            medianprops=dict(color=c3))
#for j,i in enumerate(x):
#    plt.hlines(y = result_baseline_mp_ac[j], xmin=i+25,xmax=i+175, color = 'black')
#    plt.hlines(y = result_baseline_mp_ab[j], xmin=i-175,xmax=i-25, color = 'black')
        plt.xlim(100,800*x_len+800)
        plt.title('Cross-reactivity boxplots')
        plt.xticks(x,[str(i) for i in range(x_len)])
        plt.xlabel('Genotype')
        plt.ylabel('Cross-reactivity')
        blue_patch = mpatches.Patch(color='blue', label='C protein')
        red_patch = mpatches.Patch(color='red', label='B protein')  
        green_pacth = mpatches.Patch(color = 'green',label='D protein')
        plt.legend(handles=[blue_patch, red_patch,green_patch])
#plt.hlines(y = 0.9867, xmin=100,xmax=6000, color = 'red', linestyles='--')
#plt.hlines(y = 0.9867, xmin=100,xmax=6000, color = 'blue', linestyles='--')
#plt.hlines(y = result_base_ab_v, xmin=250,xmax=6000, color = 'red', linestyles='--')
#plt.hlines(y = result_base_ac_v, xmin=250,xmax=6000, color = 'blue', linestyles='--')
#plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,],['5','10','15','20','25','30','40','50','60','70','80','90','100'])
    #plt.savefig('valid_pair.png')
        if savefile is not None:
            plt.savefig(savefile)
            return
        plt.show()
        
    elif protein3 is not None:    
        x = np.array([i for i in range(500, 500*x_len+500,500)])
        c1 = "red"
        c2 = 'blue'
        data1 = Cross_react_predict(protein1, protein2, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        data2 = Cross_react_predict(protein1, protein3, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        plt.figure()
        plt.boxplot(data1,positions=x-100,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='pink', color=c1),
            capprops=dict(color=c1),
            whiskerprops=dict(color=c1),
            flierprops=dict(color=c1, markeredgecolor=c1),
            medianprops=dict(color=c1))
        plt.boxplot(data2,positions=x+100,widths=150,notch=True, patch_artist=True,
            boxprops=dict(facecolor='cyan', color=c2),
            capprops=dict(color=c2),
            whiskerprops=dict(color=c2),
            flierprops=dict(color=c2, markeredgecolor=c2),
            medianprops=dict(color=c2))
#for j,i in enumerate(x):
#    plt.hlines(y = result_baseline_mp_ac[j], xmin=i+25,xmax=i+175, color = 'black')
#    plt.hlines(y = result_baseline_mp_ab[j], xmin=i-175,xmax=i-25, color = 'black')
        plt.xlim(100,500*x_len+500)
        plt.title('Cross-reactivity boxplots')
        plt.xticks(x,[str(i) for i in range(x_len)])
        plt.xlabel('Genotype')
        plt.ylabel('Cross-reactivity')
        blue_patch = mpatches.Patch(color='blue', label='C protein')
        red_patch = mpatches.Patch(color='red', label='B protein')  
        plt.legend(handles=[blue_patch, red_patch])
#plt.hlines(y = 0.9867, xmin=100,xmax=6000, color = 'red', linestyles='--')
#plt.hlines(y = 0.9867, xmin=100,xmax=6000, color = 'blue', linestyles='--')
#plt.hlines(y = result_base_ab_v, xmin=250,xmax=6000, color = 'red', linestyles='--')
#plt.hlines(y = result_base_ac_v, xmin=250,xmax=6000, color = 'blue', linestyles='--')
#plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,],['5','10','15','20','25','30','40','50','60','70','80','90','100'])
    #plt.savefig('valid_pair.png')
        if savefile is not None:
            plt.savefig(savefile)
            return
        plt.show()
    else:
        data1 = Cross_react_predict(protein1, protein2, genotypes, immun_predictor,cross_react_predictor,N_iterations,k,N_proc,N_individuals)
        plt.figure()
        plt.boxplot(data1)
        plt.title('Cross-reactivity boxplots')
        plt.xlabel('Genotype')
        plt.ylabel('Cross-reactivity')
        if savefile is not None:
            plt.savefig(savefile)
            return
        plt.show()


def protein_combine(*args):
    new_pool = []
    new_seq = ''
    for arg in args:
        new_seq += arg.protein_seq + '+'
        new_pool += arg.peptide_pool
    new_protein = Protein(new_seq)
    new_protein.peptide_pool = sorted(new_pool)
    return new_protein
    

    