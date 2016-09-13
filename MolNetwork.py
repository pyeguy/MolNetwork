"""A tool for generating networks based on molecular fingerprint similarity"""


__author__ = "Cameron Pye"
__email__ = "cameron.pye@gmail.com"
__license__ = "GPL"
__version__ = "0.1"
__status__ = "beta"

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

FPMETHODS = {
    'ECFC4' : lambda x:AllChem.GetMorganFingerprint(x,4,useChirality=True),
    'ECFP4' : lambda x:AllChem.GetMorganFingerprintAsBitVect(x,4,useChirality=True),
}

METRICS = {
    'DICE' : AllChem.DataStructs.DiceSimilarity,
    'Tanimoto' : AllChem.DataStructs.TanimotoSimilarity,
    'Cosine' : AllChem.DataStructs.CosineSimilarity,
    'Sokal' : AllChem.DataStructs.SokalSimilarity,
    'Tversky' : AllChem.DataStructs.TverskySimilarity,
}

BINARYFPS = {'ECFP4'}
BINONLYMETRICS = {'Cosine','Sokal','Tversky'}

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
            if mol:
                i +=1
                name = mol.GetProp('_Name')
                if name == '':
                    name = i
                yield [mol,name]


def genFps(infile,fptype):
    fpmethod = FPMETHODS[fptype]
    infile_components = infile.split('.')
    try:
        with open('{}_{}.pkl'.format(infile_components[0],fptype),'rb') as fin:
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
            fps.append([name,fpmethod(mol)])
        with open('{}_{}.pkl'.format(infile_components[0],fptype),'wb') as fout:
            pickle.dump(fps,fout)
    return fps 

def half_simtable(fps,metric):
    metric_method = METRICS[metric]
    print("Making Sim Table")
    i = 1
    scores = []
    non_sings = set()
    bar = progressbar.ProgressBar()
    for name1,fp1 in bar(fps):
        for name2,fp2 in fps[i:]:
            # if j > i :
            score = metric_method(fp1,fp2)
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

def directed_simtable(fps,weights):
    metric_method = METRICS['Tversky']
    print("Making Full Directed Sim Table")
    i = 1
    j = 1
    scores = []
    non_sings = set()
    bar = progressbar.ProgressBar()
    for name1,fp1 in bar(fps):
        for name2,fp2 in fps:
            if i != j:
                score = metric_method(fp1,fp2)
                if score >= args.thresh:
                    non_sings.add(name1)
                    non_sings.add(name2)
                    scores.append([name1,name2,score])
            j+=1
        i+=1
    names = set([x[0] for x in fps])
    sings = names - names.intersection(non_sings)
    print('\n','-'*80,sep='')
    print("Percent clustered : {0:.2f}%".format(((len(names)-len(sings))/len(names)*100)))
    print("Number of Singletons : {}".format(len(sings)))
    if args.singletons:
        print("including singletons in edgelist as self loops")
        scores.extend([[sing,sing,1.0] for sing in sings])
    return scores

def sample_scores(fps,thresh,samplesize,metric):
    metric_method = METRICS[metric]
    if samplesize > len(fps):
        raise Exception('Sample Size cannot exceede number of molecules')
    scores = []
    non_sings = set()
    pbar = progressbar.ProgressBar()
    for name1,fp1 in pbar(fps):
        for name2,fp2 in sample(fps,samplesize):
            if fp1 != fp2:
                score = metric_method(fp1,fp2)
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

def negSAR(full_filename,hits_fps,metric):
    metric_method = METRICS[metric]
    # full_filename_components = full_filename.split('.')
    neg_sar_scores = []
    neg_sar_names = set()
    # full_fps = {name:fp for name,fp in genFps(full_mollist)}
    print('Getting Fingerprints of the full library... ')
    full_fps = genFps(full_filename,args.fptype)
    pbar = progressbar.ProgressBar()
    print('Making Neg SAR comparisons')
    for hit_name,hit_fp in pbar(hits_fps):
        for neg_name,neg_fp in full_fps:
            if hit_name != neg_name:
                score = metric_method(hit_fp,neg_fp)
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
            * ...
        '''))
    parser.add_argument('infile',help='input file; smiles file for now')
    parser.add_argument('outfile',help='output edgelist file')
    parser.add_argument('--thresh','-t',default=0.8,help='threshold for edge cutoff' , type=float)
    parser.add_argument('--samplesize','-s',help='number of compounds to sample when generating the sim table (for large # of compounds)',type=int)
    parser.add_argument('--singletons',action='store_true',help='include self-loop for all singletons')
    parser.add_argument('--negsar',help='File with all compounds to be clustered w/hits. will output a second file with hit attaribute')
    parser.add_argument('--fptype',help='the type of fp to be used (these are the Morgan Algo. equivalents)',choices=['ECFC4','ECFP4'],default='ECFC4')
    mutex_group1 = parser.add_mutually_exclusive_group()
    mutex_group1.add_argument('--metric',help='Distance/Simiarity Metric to use',choices=['DICE','Tanimoto','Cosine','Sokal'],default='DICE')
    mutex_group1.add_argument('--directed',help='Use Tversky scoring to generate directed graph. include when using this option <alpha> <beta>',type=list)
    parser.add_argument('--nx', help='Output a .gpkl NetworkX graph pickle file',action='store_true')
    parser.add_argument('--gexf',help='Output a .gexf file which can be read by Gephi natively or Cytoscape with the gexf-app plugin',action='store_true')
    args = parser.parse_args()

    infile_components = args.infile.split('.')
    outfile_components = args.outfile.split('.')


    if (args.metric in BINONLYMETRICS and args.fptype not in BINARYFPS) or (args.directed and args.fptype not in BINARYFPS):
        parser.print_help()
        print()
        if args.directed:
            raise Exception("Unspported Metric and FP combination : for directed use binary Fingerprint ie --fptype ECFP4")
        else: raise Exception("Unspported Metric and FP combination : <{}> & <{}>".format(args.metric,args.fptype))


    fps = genFps(args.infile,args.fptype)
    
    if args.samplesize:
        scores = sample_scores(fps,args.thresh,args.samplesize,args.metric)
    else:
        scores = half_simtable(fps,args.metric)

    if args.negsar:
        neg_sar_scores, neg_sar_names = negSAR(args.negsar,fps,args.metric)
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
    if args.nx or args.gexf:
        import networkx
        if args.directed:
            g = networkx.DiGraph()
        else:
            g = networkx.Graph()
        ebunch = [(s,t,{'weight':w}) for s,t,w in scores]
        g.add_edges_from(ebunch)
        if args.nx:
            networkx.write_gpickle(g,'{}.gpkl'.format(outfile_components[0]))
        if args.gexf:
            networkx.write_gexf(g,'{}.gexf'.format(outfile_components[0]))