import yasara, sys, os, time, multiprocessing
import numpy as np
import pandas as pd

# USAGE (after loading the corresponding conda environment): python E91D.py system.sce

def yasaraInitialization():
    yasara.info.mode = "txt"
    yasara.Processors(8,0)

def loadMDSimulation(scene):
    # get sim file names and associated snapshot ids
    simname, simnum = [], []
    for elem in os.listdir("../../../traj-ABCD/PDBs-wat-wtABCD300/"):
        if ".pdb" in elem:
            simname.append(elem)
    # metrics to be analized
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    for frame in simname:
        pool.apply_async(surfingSnapshots, args=(scene,frame,))
    pool.close()
    pool.join()

def surfingSnapshots(scene,frame):
    # load sim files, remove cell and water, save them in yob
    yasaraInitialization()
    yasara.LoadSce(scene)
    #yasara.Sim(control = "Pause")
    #yasara.LoadSim(frame)
    yasara.LoadPDB("../../../traj-ABCD/PDBs-wat-wtABCD300/" + frame)
    simnum = int(frame[-8:-4])
    #yasara.DelObj("2 3")
    #delete the obj corresponding to the loaded sce
    yasara.DelObj("1")
    # delete xtal waters in obj 1
    #yasara.DelRes("water cip cim")
    # calculate various types of features
    hbonds,hbondD94E,hbondG93,hbondL95,hbondH121,hbondR122,hbondsol,salt,saltR122 = calculateInteraction() 
    yasara.Exit()
    df = pd.DataFrame({"Frame":[simnum], "Hbonds":hbonds,"HbondD94E":hbondD94E,"HbondG93":hbondG93,"HbondL95":hbondL95,"HbondH121":hbondH121, \
                             "HbondR122":hbondR122,"Hbonds solv":hbondsol, \
                             "Salt":salt,"SaltR122":saltR122}, \
                          columns=["Frame", "Hbonds","HbondD94E", "HbondG93", "HbondL95", "HbondH121", \
                          "HbondR122", "Hbonds solv", "Salt", "SaltR122"]) 
    df.to_csv("tmpAdvancedStructuralAnalysisMultiprocess" + str(simnum) + ".csv", index=False)
    print(df)

def calculateInteraction():
    # protein hbonds
    #I divide by 4 because it is a tetramer
    hbonds = [len(yasara.ListHBoAtom("protein and res 89 375 661 947","protein",Min=0,results=1))/4]
    hbondD94E = [len(yasara.ListHBoAtom("protein and res 89 375 661 947","protein and res 92 378 664 950",Min=0,results=1))/4]
    hbondG93 = [len(yasara.ListHBoAtom("protein and res 89 375 661 947","protein and res 91 377 663 949",Min=0,results=1))/4]
    hbondL95 = [len(yasara.ListHBoAtom("protein and res 89 375 661 947","protein and res 93 379 665 951",Min=0,results=1))/4]
    hbondH121 = [len(yasara.ListHBoAtom("protein and res 89 375 661 947","protein and res 119 405 691 977",Min=0,results=1))/4]
    hbondR122 = [len(yasara.ListHBoAtom("protein and res 89 375 661 947","protein and res 120 406 692 978",Min=0,results=1))/4]

    # protein-solvent hbonds
    hbondsol = [len(yasara.ListHBoAtom("protein and res 89 375 661 947","res WAT",Min=0,results=1))/4]    

    # protein-Na+ interaccions (better atom here)
    #sodium = [len(yasara.ListConAtom("res Na+","element oxygen protein and res 118",cutoff=3.5,exclude=4,occluded="no",sort="No"))]  

    # protein salt bridges
    salt = [len(yasara.ListIntRes("sidechain and res 89 375 661 947","sidechain",Type="Ionic",cutoff=5.0,exclude=4,occluded="Yes",sort="No"))//4]
    saltR122 = [len(yasara.ListIntRes("sidechain and res 89 375 661 947","sidechain and res 120 406 692 978",Type="Ionic",cutoff=5.0,exclude=4,occluded="Yes",sort="No"))//4]

    # protein hydrophobic contacts counting atoms in sc or only one count per sidechain
    #phobiclist = [len(yasara.ListIntRes("sidechain res Ile Leu Val Ala Phe Pro Trp Tyr Met","sidechain res Ile Leu Val Ala Phe Pro Trp Tyr Met",Type="Hydrophobic",cutoff=5.0,exclude=4,occluded="No",sort="No"))]
    #protein pipi at NON-canonical interfaces
    #pipiAB = [len(yasara.ListIntRes("sidechain res Phe Tyr Trp and res 1-283","sidechain res Phe Tyr Trp and res 287-569",Type="PiPi",cutoff=5.0,exclude=4,occluded="Yes",sort="No"))]


    return hbonds,hbondD94E,hbondG93,hbondL95,hbondH121,hbondR122,hbondsol,salt,saltR122

def concatenateSnapshotsAnalysis():
    dataframes = []
    for elem in os.listdir("./"):
        if "tmpAdvancedStructuralAnalysisMultiprocess" in elem:
            dataframes.append(pd.read_csv(elem))
    df = pd.concat(dataframes)
    df.sort_values(["Frame"], ascending=True, inplace=True)
    df.to_csv("E91D-wtABCD300.csv", index=False)
    #os.system("mkdir tmp")
    #os.system("mv tmpAdvancedStructuralAnalysisMultiprocess* ./tmp")
    os.system("rm tmpAdvancedStructuralAnalysisMultiprocess*")

#### MAIN ####

start = time.time()
loadMDSimulation(sys.argv[1])
concatenateSnapshotsAnalysis()
print("time (min): ", (time.time() - start)/60.)
yasara.Exit()
