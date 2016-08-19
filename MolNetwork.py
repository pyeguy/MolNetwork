from rdkit import Chem
from rdkit.Chem import AllChem
import progressbar
from random import sample
import pickle
import argparse
import textwrap
# import itertools
# from sys import argv
# script,infile,outfile = argv

#list comp versions of the for loops.
#fps = [[AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi),4,useChirality=True),name] for smi,name in smilist]
#scores = [[name1,name2,Chem.DataStructs.DiceSimilarity(fp1,fp2)] for name1,fp1 in fps for name2,fp2, in fps]



def getNummols(infile):
    infile_components = infile.split('.')
    nummols = 0
    if infile_components[1] == 'smi':
        with open(infile) as fin:
            for line in fin:
                nummols +=1
    elif infile_components[1] in {'sdf','sd'}:
        suppl = Chem.SDMolSupplier(infile)
        nummols = len(suppl)
    return nummols

def genMollist(infile):
    infile_components = infile.split('.')
    if infile_components[1] == 'smi':
        with open(infile) as fin:
            for line in fin: 
                line = line.strip()
                if line:
                    line = line.split('\t')
                    yield [Chem.MolFromSmiles(line[0]),line[1]]
    elif infile_components[1] in {'sdf','sd'}:
        suppl = Chem.SDMolSupplier(infile)
        i = 0
        for mol in suppl:
            i +=1
            name = mol.GetProp('_Name')
            if name == '':
                name = i
            yield [mol,name]


def genFps(infile):
    infile_components = infile.split('.')
    try:
        with open(infile_components[0]+'.pkl','rb') as fin:
            print('Loading old FP pkl file')
            fps = pickle.load(fin)
            print('Done Loading FPs')
    except:
        nummols = getNummols(infile)
        molgen = genMollist(infile)
        print("Making FPs")
        fps = []
        pbar = progressbar.ProgressBar(max_value=nummols)
        for mol,name in pbar(molgen):
            fps.append([name,AllChem.GetMorganFingerprint(mol,4,useChirality=True)])
        with open(infile_components[0]+'.pkl','wb') as fout:
            pickle.dump(fps,fout)
    return fps 

def full_simtable(fps):
    print("Making Sim Table")
    i = 1
    scores = []
    non_sings = set()
    bar = progressbar.ProgressBar()
    for name1,fp1 in bar(fps):
        for name2,fp2 in fps[i:]:
            # if j > i :
            score = Chem.DataStructs.DiceSimilarity(fp1,fp2)
            if score >= args.thresh:
                non_sings.add(name1)
                non_sings.add(name2)
                scores.append([name1,name2,score])
        i +=1
    names = set([x[0] for x in fps])
    sings = names - names.intersection(non_sings)
    print('\n','-'*80,sep='')
    print("Percent clustered : {0:.2f}%".format(((len(names)-len(sings))/len(names)*100)))
    print("Number of Singletons : {}".format(len(sings)))
    if args.singletons:
        print("including singletons in edgelist as self loops")
        scores.extend([[sing,sing,1.0] for sing in sings])
    return scores

def sample_scores(fps,thresh,samplesize):
    if samplesize > len(fps):
        raise Exception('Sample Size cannot exceede number of molecules')
    scores = []
    non_sings = set()
    pbar = progressbar.ProgressBar()
    for name1,fp1 in pbar(fps):
        for name2,fp2 in sample(fps,samplesize):
            if fp1 != fp2:
                score = Chem.DataStructs.DiceSimilarity(fp1,fp2)
                if score >= thresh:
                    non_sings.add(name1)
                    non_sings.add(name2)
                    scores.append([name1,name2,score])
    names = set([x[0] for x in fps])
    sings = names - names.intersection(non_sings)
    print('\n','-'*80,sep='')
    print("Percent clustered : {0:.2f}%".format(((len(names)-len(sings))/len(names)*100)))
    print("Number of Singletons : {}".format(len(sings)))
    if args.singletons:
        print("including singletons in edgelist as self loops")
        scores.extend([[sing,sing,1.0] for sing in sings])  
    return scores

def negSAR(full_filename,hits_fps):
    # full_filename_components = full_filename.split('.')
    neg_sar_scores = []
    neg_sar_names = set()
    # full_fps = {name:fp for name,fp in genFps(full_mollist)}
    print('Getting Fingerprints of the full library... ')
    full_fps = genFps(full_filename)
    pbar = progressbar.ProgressBar()
    print('Making Neg SAR comparisons')
    for hit_name,hit_fp in pbar(hits_fps):
        for neg_name,neg_fp in full_fps:
            if hit_name != neg_name:
                score = Chem.DataStructs.DiceSimilarity(hit_fp,neg_fp)
                if score >= args.thresh:
                    neg_sar_scores.append([hit_name,neg_name,score])
                    neg_sar_names.add(neg_name)
    return neg_sar_scores,neg_sar_names



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class= argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Generate a network edge list from a molecule file
        -------------------------------------------------
        Uses Morgan Feature Counts and Dice Similarity Scores 
        
        Input file can be .smi file or .sd file
            * smi file should be of the format : smiles<tab>name,

        Output file(s) Will be generated using the output file name
            * A pkl file will be generated storing the FPs for future use
            * Specified output name will be an edgelist file (source<tab>targe<tab>score,)
            * If --negsar option is used then a hit map will also be generated
                * 1 = hit , 0 = a compound that clustered with a hit but did not hit
        
        ** Experimental **
        If a very large input file is used (ie >10k molecules) generating the full 
        pairwise scores should not be attempted. Instead try using the --sample <n> argument. 
        This will sample n compounds from the whole to attempt to psuedo cluster the compounds

        Future Versions Will Include:
            * multiprocessing
            * richer input output options
                * networkx / gephi output
            * ...
        '''))
    parser.add_argument('infile',help='input file; smiles file for now')
    parser.add_argument('outfile',help='output edgelist file')
    parser.add_argument('--thresh','-t',default=0.8,help='threshold for edge cutoff' , type=float)
    parser.add_argument('--samplesize','-s',help='number of compounds to sample when generating the sim table (for large # of compounds)',type=int)
    parser.add_argument('--singletons',action='store_true',help='include self-loop for all singletons')
    parser.add_argument('--negsar',help='File with all compounds to be clustered w/hits. will output a second file with hit attaribute')
    args = parser.parse_args()

    infile_components = args.infile.split('.')
    outfile_components = args.outfile.split('.')

    fps = genFps(args.infile)
    
    if args.samplesize:
        scores = sample_scores(fps,args.thresh,args.samplesize)
    else:
        scores = full_simtable(fps)

    if args.negsar:
        neg_sar_scores, neg_sar_names = negSAR(args.negsar,fps)
        print('Neg SAR compounds: {}'.format(len(neg_sar_names)))
        scores.extend(neg_sar_scores)
        hitmap = [[neg,0] for neg in neg_sar_names]
        hitmap.extend([hit,1] for hit in [x[0] for x in fps])
        with open('{}_hitMap.txt'.format(outfile_components[0]),'w') as fout:
            print('Name\tHitBin',file=fout)
            for line in hitmap:
                print('{}\t{}'.format(*line),file=fout)

    with open(args.outfile,'w') as fout:
        print('Source\tTarget\tScore',file=fout)
        for intline in scores:
            print('{}\t{}\t{}'.format(*intline),file=fout)
