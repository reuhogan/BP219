from pyrosetta import *
init()
from pyrosetta.toolbox import *
from pyrosetta import PyMOLMover
import sys
pmm=PyMOLMover()
import pandas as pd
from pyrosetta.rosetta.core.scoring import CA_rmsd

def mutateresis(pose):
    fa_sfxn=create_score_function("ref2015")
    fastrelax=pyrosetta.rosetta.protocols.relax.FastRelax()
    fastrelax.set_scorefxn(fa_sfxn)

    resnums_w_mutres={} #define a dictionary with all PDB positions and mutations to be made
    resnums_w_mutres[484]=['E','Q']
    resnums_w_mutres[452]=['L','R']
    resnums_w_mutres[478]=['T','K']

    mutname,score,rmsd=[],[],[]
    for respos,resids in zip(resnums_w_mutres.keys(),resnums_w_mutres.values()): #iterates through dictionary terms
        origres, mutateto=resids[0], resids[1]
        position=int(pose.pdb_info().pdb2pose("E",respos)) 
        mutate_residue(temppose, position, mutateto) #mutate the listed residue and chain E (taken from pdb)
        fastrelax.apply(temppose)
        pmm.pymol_name("%s%s%s_pose" % (origres,position,mutateto))
        pmm.apply(temppose)
        temppose.dump_pdb("bp219/%s%s%s_pose.pdb" % (origres,position,mutateto)) #output the mutated residue to files
        fascore=fa_sfxn.score(temppose)
        score.append(fascore)
        mutname.append("%s%s%s" % (origres,position,mutateto))
        rmsd.append(CA_rmsd(temppose, pose))
    return mutname,score,rmsd

def preprocessing(origposepath):
    pose = pose_from_pdb(origposepath) #initialize the pdb into a pose object
    fa_sfxn=create_score_function("ref2015") #initialize the rosetta score function - using a full atom scoring method 
    fastrelax=pyrosetta.rosetta.protocols.relax.FastRelax() #instance of a energy minimization function
    fastrelax.set_scorefxn(fa_sfxn)
    fastrelax.apply(pose)
    origposescore=fa_sfxn.score(pose)
    return pose,origposescore #return pose object and its rosetta energy score

def main():
    origposepath=sys.argv[1] #accepts the path to the starting spike protein 
    pose,origposescore=preprocessing(origposepath)
    mutname,score,rmsd=mutateresis(pose)
    mutname.append("orignal_structure") #adds on data from starting structure to the lists that will comprise the dataframe
    score.append(origposescore) #Rosetta score of the starting structure
    rmsd.append(0) #RMSD will be 0 
    outputinfo=pd.DataFrame()#table to keep track of all output's scores
    for i,j in zip(["Mutations","Rosetta_Scores","RMSD"],[mutname,score,rmsd]):
        outputinfo[i]=j
    outputinfo.to_csv("OutputInfo.csv") #output the data to a CSV for further processing

if __name__ == '__main__':
    main()    


