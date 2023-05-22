#from ._version import __version__
#from . import abcrpv
from . import rpv_definitions as rdef
from . import rpv_misc as rmisc

import numpy as np
import pandas as pd  
from IPython.display import display 
import subprocess
import itertools
import warnings
import os 
from pathlib import Path

abcrpv_package_path = Path(__file__).parent.absolute()


if not os.path.isfile(os.path.join(abcrpv_package_path,"input/table_notsup.csv")):
    raise FileNotFoundError("Input transition tables {x} does not exist, check paths or download from git repository".format(x=os.path.join(abcrpv_package_path,"input/table_notsup.csv")))

if not os.path.isfile(os.path.join(abcrpv_package_path,"input/table_sup.csv")):
    raise FileNotFoundError("Input transition tables {x} does not exist, check paths or download from git repository".format(x=os.path.join(abcrpv_package_path,"input/table_sup.csv")))

if not os.path.isfile(os.path.join(abcrpv_package_path,"input/table_strsup.csv")):
    raise FileNotFoundError("Input transition tables {x} does not exist, check paths or download from git repository".format(x=os.path.join(abcrpv_package_path,"input/table_strsup.csv")))

NOTSUP_TABLE   = pd.read_csv(os.path.join(abcrpv_package_path,"input/table_notsup.csv"))
SUP_TABLE      = pd.read_csv(os.path.join(abcrpv_package_path,"input/table_sup.csv"))
STRSUP_TABLE   = pd.read_csv(os.path.join(abcrpv_package_path,"input/table_strsup.csv"))
SUPDICT= {
 "notsup"             : NOTSUP_TABLE,
 "notsuppressed"      : NOTSUP_TABLE,
 "sup"                : SUP_TABLE,
 "suppressed"         : SUP_TABLE,
 "strsup"             : STRSUP_TABLE,    
 "stronglysup"        : STRSUP_TABLE,    
 "strsuppressed"      : STRSUP_TABLE,    
 "stronglysuppressed" : STRSUP_TABLE,    
}

def transition_sig(mother,daughter,sup,opt=False):
    """
    Return signatures from the transitions of mother sparticle to daughter sparticles
    @TODO ref to code
    
    Inputs
    ------
    -mother    (str): mother sparticle   (should be in sparticles-format)      
    -daugther  (str): daughter sparticle (should be in sparticles-format)            
    -sup       (str): suppression level  ("sup", "notsup", "strsup")
                      the suppression level of transitions is based on studies in arxiv:1205.0557. It has been modified and should refer to:
                      - input/table_notsup.dat
                      - input/table_strsup.dat
                      - input/table_sup.dat
    -opt            : print message if cant find it(default=False)
    """
    trans_pd = SUPDICT[sup.lower().replace(" ","")]
    trans_pd[(trans_pd["Mother"] == rmisc.check_format(mother)) & (trans_pd["Daughter"] == rmisc.check_format(daughter))]
    transition = trans_pd[(trans_pd["Mother"] == rmisc.check_format(mother)) & (trans_pd["Daughter"] == rmisc.check_format(daughter)) ]["Transitions"].values
    if len(transition) == 0:
        if opt == True:
            print("\n"+(rmisc.check_format(mother)+" to "+rmisc.check_format(daughter)+" is not "+sup))
        return None
    assert len(transition) == 1, ("Should be one occurance of the transition in table_{sup}.csv".format(sup=sup), transition)
    return transition[0].split()

def set_elements_simplified(instr):
    """
    In terms of signatures:
        Get all possible signatures from signatures with simplified objects (J,3,L)
    or in terms of sets:
        Get all elements from signature set defined via simplified objects (J,3,L)

    @TODO SF bracket treament

    Inputs
    ------
    -instr  (str): input signature in one-char syntax  

    Example
    ------
    JJjtLTEE => ['tjjjTlE', 'tjjjTTE', 'ttjjTlE', 'ttjjTTE', 'tbjjTlE', 'tbjjTTE', 'ttjjTlE', 'ttjjTTE', 'tttjTlE', 'tttjTTE', 'ttbjTlE', 'ttbjTTE', 'tbjjTlE', 'tbjjTTE', 'ttbjTlE', 'ttbjTTE', 'tbbjTlE', 'tbbjTTE']
    """
    assert rmisc.check_signature_format(instr) == True, "Check Signature Format"
    strlist = list(instr)
    expand_list = []
    for i in strlist:
        if i == "J":
            expand_list.append(["j","t","b"])
        elif i == "3":
            expand_list.append(["t","b"])
        elif i == "L":
            expand_list.append(["l","T"])
        else:
            expand_list.append(i)
    return [rmisc.signature_ordering("".join(x)) for x in list(itertools.product(*expand_list))]

def get_all_superset(instr,optJ,opt3,optL):
    """
    In terms of signatures:
        Get all possible simplified signatures that can lead to input signature 
    or in terms of sets:
        Get all possible superset (defined via simplified objects (J,3,L)) that contain instr
    @TODO SF bracket treament

    Inputs
    ------
    -instr  (str) : input signature in one-char syntax  
    -optX   (bool): Flag to look up simplified signatures with X

    Example
    ------
    Input = tjTlE
      optJ   opt3    optL
      True   True    True  ['JJLLE', 'JJLlE', 'JJLTE', 'JJTlE', 'JjLLE', 'JjLlE', 'JjLTE', 'JjTlE', 'J3LLE', 'J3LlE', 'J3LTE', 'J3TlE', '3jLLE', '3jLlE', '3jLTE', '3jTlE', 'JtLLE', 'JtLlE', 'JtLTE', 'JtTlE', 'tjLLE', 'tjLlE', 'tjLTE', 'tjTlE']
      True   False   False ['JJTlE', 'JjTlE', 'JtTlE', 'tjTlE']
      False  True    False ['3jTlE', 'tjTlE']
      False  False   True  ['tjLLE', 'tjLlE', 'tjLTE', 'tjTlE']
      True   True    False ['JJTlE', 'JjTlE', 'J3TlE', '3jTlE', 'JtTlE', 'tjTlE']
      True   False   True  ['JJLLE', 'JJLlE', 'JJLTE', 'JJTlE', 'JjLLE', 'JjLlE', 'JjLTE', 'JjTlE', 'JtLLE', 'JtLlE', 'JtLTE', 'JtTlE', 'tjLLE', 'tjLlE', 'tjLTE', 'tjTlE']
      False  True    True  ['3jLLE', '3jLlE', '3jLTE', '3jTlE', 'tjLLE', 'tjLlE', 'tjLTE', 'tjTlE']
    """
    assert rmisc.check_signature_format(instr) == True, "Check Signature Format"
    strlist = list(instr)
    expand_list = []
    for i in strlist:
        if i == "j" and optJ==True:
            expand_list.append(["J",i])
        elif i == "t" or i == "b":
            if optJ==True and opt3==True:
                expand_list.append(["J","3",i])
            elif optJ==True:
                expand_list.append(["J",i])
            elif opt3==True:
                expand_list.append(["3",i])
            else:
                expand_list.append(i)    
        elif i == "l" and optL==True :
            expand_list.append(["L",i])
        elif i == "T" and optL==True:
            expand_list.append(["L",i])
        else:
            expand_list.append(i)
    return [rmisc.signature_ordering("".join(x)) for x in list(itertools.product(*expand_list))]

def get_subsets(insiglist):
    """
    In terms of signatures:
        Given a list of signatures, get all possible simplified signatures that contains them
    or in terms of sets:
        Given the unniverse set, get all possible subsets

    Inputs
    ------
    -insiglist    (array/list of str) : array/list of signatures in one-char syntax
    
    Output:
    ------
    An array of [[subset_1_name, (elements of subset_1)],[subset_2_name, (elements of subset_2)],...]
    """
    if type(insiglist) == str:
        print("insiglist is a single string")
        insiglist=[insiglist]

    #input is list of signatures that should be the universe set
    insiglist = set(insiglist)

    #make sure everything is ordered
    insiglist = (list(map(rmisc.signature_ordering,insiglist)))

    #get all possible elements of universe set
    insiglist = (list(map(set_elements_simplified,insiglist)))
    insiglist = list(itertools.chain(*insiglist))
    insiglist = set(insiglist)

    #Group via number of objects (J...L...MET)
    subsets = []
    sig_J_L =  set(list(map(rmisc.to_L,set(list(map(rmisc.to_J,insiglist))))))
    for i in sig_J_L:
        x = [s for s in insiglist if ((rmisc.to_L((rmisc.to_J(s))))) == i]
        cJ = rmisc.check_merge(x,"J")
        c3 = rmisc.check_merge(x,"3")
        cL = rmisc.check_merge(x,"L")
        for j in x:
            for k in get_all_superset(j,cJ,c3,cL):
                if set(set_elements_simplified(k)).issubset(insiglist):
                    if [k,set(set_elements_simplified(k))] not in subsets:
                        subsets.append([k,set(set_elements_simplified(k))])
    # Return an array of [subset_name, (elements of subset)]
    return subsets

def minimal_sets_simplified_signatures(insiglist):
    """
    Get the minimal subsets required to cover all signatures in insiglist
    
    Inputs
    ------
    -insiglist    : list of signatures
    
    Output:
    ------
    -output     : 
    -noutput    :  
    """    

    if type(insiglist) == str:
        print("insiglist is stirng")
        insiglist = [insiglist]

    #get all possible subsets that can provide the signatures (get_subsets)
    outtemp = get_subsets(insiglist)

    #from all possible subsets, find minimal sets required to get all siagntures (setCover))
    tempcover = rmisc.setCover(list(np.transpose(outtemp)[1]))

    output = []
    noutput = []
    for j in tempcover:
        noutput.append(len(j))
        output.append([x[0] for x in outtemp if x[1] == j][0])
    return output,noutput

"""
def generate_transition_df(save_csv=False):
    #Get transition tables
    dat = []
    for i in rdef.SPARTICLES:
        print(".",end="")
        for j in rdef.SPARTICLES:
            if i == j :
                continue
            
            #1vertex transition (i-j)   
            if transition_sig(i,j,"notsup"):
                dat.append([i,j,1,"-",rmisc.signature_ordering(transition_sig(i,j,"notsup")[0]),i+" -- ("+transition_sig(i,j,"notsup")[0]+") -- "+j])

            #n-vertex chain    
            else:
                suptemp = []

                #2-vertex transition (i-k-j)    
                for k in rdef.SPARTICLES:
                    if k != i and k != j and transition_sig(i,k,"notsup") and transition_sig(k,j,"notsup"):
                        dat.append([i,j,2,k,rmisc.signature_ordering(transition_sig(i,k,"notsup")[0]+transition_sig(k,j,"notsup")[0]),i+" -- ("+transition_sig(i,k,"notsup")[0]+") -- "+k+" -- ("+transition_sig(k,j,"notsup")[0]+") -- "+j])
                        suptemp.append([k,transition_sig(i,k,"notsup"),transition_sig(k,j,"notsup")])

                #3-vertex transition (i-k-l-j)    
                if len(suptemp) == 0:
                    for k in rdef.SPARTICLES:
                        for l in rdef.SPARTICLES:
                            if i != k and k != l and l != j and  i != l and k != j and transition_sig(i,k,"notsup") and transition_sig(k,l,"notsup") and transition_sig(l,j,"notsup"):
                                dat.append([i,j,3,k+" "+l,rmisc.signature_ordering(transition_sig(i,k,"notsup")[0]+transition_sig(k,l,"notsup")[0]+transition_sig(l,j,"notsup")[0]),i+" -- ("+transition_sig(i,k,"notsup")[0]+") -- "+k+" -- ("+transition_sig(k,l,"notsup")[0]+") -- "+l+" -- ("+transition_sig(l,j,"notsup")[0]+") -- "+j])
                                suptemp.append([k,transition_sig(i,k,"notsup"),transition_sig(k,j,"notsup")])

                if len(suptemp) == 0:    
                    print(i,j,": More than a 3 vertex decay")

            #some checks for one obj sup,strsup decay
            #if i == "H^0" or j == "H^0" or i == "H^+" or j == "H^+":
            #    if transition_sig(i,j,"H"):
            #        if len(transition_sig(i,j,"H")[0]) == 1:
            #            dat.append([i,j,1,"-",transition_sig(i,j,"H")[0],i+" -- ("+transition_sig(i,j,"H")[0]+") -- "+j])
    
    index_name = ['Mother', 'Daughter', 'Nvertex', 'Intermediate', 'Signatures', 'Chain']
    transition_df= pd.DataFrame(columns=index_name)
    entry = pd.DataFrame(dat, columns=index_name)
    transition_df = transition_df.append(entry)
    #if intermediate_csv == True:
    #    transition_df.to_csv('CSV/transition_df.csv',index=False)
    
    transition_df["Signatures (Easy-Read)"] = transition_df["Signatures"].apply(rmisc.easy_read)
    #if intermediate_csv == True:
    #    transition_df.to_csv('CSV/transition_df_EasyRead.csv',index=False)
    
    transition_df["All Possible Signature"] = ""
    for i in rdef.SPARTICLES:
        for j in rdef.SPARTICLES:
            transition_df.loc[(transition_df['Mother'] == i) & (transition_df['Daughter'] == j), ['All Possible Signature']] = set(list(map(rmisc.signature_ordering, np.array(transition_df.loc[(transition_df['Mother'] == i) & (transition_df['Daughter'] == j) ]['Signatures']))))
    #if intermediate_csv == True:
    #    transition_df.to_csv('CSV/transition_df_All_Possible_Signature.csv',index=False)

    transition_df["Minimal Set"] = transition_df["All Possible Signature"].apply(minimal_sets_simplified_signatures).str[0]
    transition_df["Minimal Set (number of elements)"] = transition_df["All Possible Signature"].apply(minimal_sets_simplified_signatures).str[1]

    #display(transition_df)
    #if intermediate_csv == True:
    #    transition_df.to_csv('CSV/transition_df_minimal_set.csv',index=False)
    
    transition_df["Minimal Set (same flavour)"] = transition_df["Minimal Set"].apply(rmisc.same_flavour)

    if save_csv == True:
        #transition_df.to_csv('CSV/transition_df_minimal_set_wf.csv',index=False)
        transition_df.to_csv(os.path.join(abcrpv_package_path,"data/main.csv"),index=False)
    print()
    return transition_df
"""

def generate_transition_df(nv=1,save_csv=False):
    #Get transition tables
    dat = []
    for i in rdef.SPARTICLES:
        print(".",end="")
        for j in rdef.SPARTICLES:
            # Ignore transitions to oneself
            if i == j :
                continue

            suptemp = False
            n1_transition = transition_sig(i,j,"notsup")

            #1vertex transition (i-j)   
            #(i,j,n1_transition,not n1_transition,nv > 1)
            if n1_transition:
                dat.append([i,j,1,"-",rmisc.signature_ordering(n1_transition[0]),i+" -- ("+n1_transition[0]+") -- "+j])

            #2-vertex transition (i-k-j)    
            if (not n1_transition) or (nv > 1) :
                for k in rdef.SPARTICLES:
                    if k != i and k != j and transition_sig(i,k,"notsup") and transition_sig(k,j,"notsup"):
                        dat.append([i,j,2,k,
                                    rmisc.signature_ordering(transition_sig(i,k,"notsup")[0]+transition_sig(k,j,"notsup")[0]),
                                    i+" -- ("+transition_sig(i,k,"notsup")[0]+") -- "+k+" -- ("+transition_sig(k,j,"notsup")[0]+") -- "+j])
                        suptemp = True
                        #suptemp.append([k,transition_sig(i,k,"notsup"),transition_sig(k,j,"notsup")])

            #3-vertex transition (i-k-l-j)
            if ((not n1_transition) and (not suptemp)) or (nv > 2) :
                    for k in rdef.SPARTICLES:
                        for l in rdef.SPARTICLES:
                            if i != k and k != l and l != j and  i != l and k != j and transition_sig(i,k,"notsup") and transition_sig(k,l,"notsup") and transition_sig(l,j,"notsup"):
                                dat.append([i,j,3,k+" "+l,rmisc.signature_ordering(transition_sig(i,k,"notsup")[0]+transition_sig(k,l,"notsup")[0]+transition_sig(l,j,"notsup")[0]),i+" -- ("+transition_sig(i,k,"notsup")[0]+") -- "+k+" -- ("+transition_sig(k,l,"notsup")[0]+") -- "+l+" -- ("+transition_sig(l,j,"notsup")[0]+") -- "+j])
                                suptemp = True
                                #suptemp.append([k,transition_sig(i,k,"notsup"),transition_sig(k,j,"notsup")])

            if not suptemp and (nv > 1) and (not n1_transition):    
                print(i,j,": More than a 3 vertex decay")
                ##@TODO: Rewrite function with a wrapping function for nvertex: 

            #some checks for one obj sup,strsup decay
            #if i == "H^0" or j == "H^0" or i == "H^+" or j == "H^+":
            #    if transition_sig(i,j,"H"):
            #        if len(transition_sig(i,j,"H")[0]) == 1:
            #            dat.append([i,j,1,"-",transition_sig(i,j,"H")[0],i+" -- ("+transition_sig(i,j,"H")[0]+") -- "+j])

    index_name = ['Mother', 'Daughter', 'Nvertex', 'Intermediate', 'Signatures', 'Chain']
    transition_df= pd.DataFrame(columns=index_name)
    entry = pd.DataFrame(dat, columns=index_name)
    transition_df = transition_df.append(entry)
    #if intermediate_csv == True:
    #    transition_df.to_csv('CSV/transition_df.csv',index=False)
    
    transition_df["Signatures (Easy-Read)"] = transition_df["Signatures"].apply(rmisc.easy_read)
    #if intermediate_csv == True:
    #    transition_df.to_csv('CSV/transition_df_EasyRead.csv',index=False)
    
    transition_df["All Possible Signature"] = ""
    for i in rdef.SPARTICLES:
        for j in rdef.SPARTICLES:
            transition_df.loc[(transition_df['Mother'] == i) & (transition_df['Daughter'] == j), ['All Possible Signature']] = set(list(map(rmisc.signature_ordering, np.array(transition_df.loc[(transition_df['Mother'] == i) & (transition_df['Daughter'] == j) ]['Signatures']))))
    #if intermediate_csv == True:
    #    transition_df.to_csv('CSV/transition_df_All_Possible_Signature.csv',index=False)

    transition_df["Minimal Set"] = transition_df["All Possible Signature"].apply(minimal_sets_simplified_signatures).str[0]
    transition_df["Minimal Set (number of elements)"] = transition_df["All Possible Signature"].apply(minimal_sets_simplified_signatures).str[1]

    #display(transition_df)
    #if intermediate_csv == True:
    #    transition_df.to_csv('CSV/transition_df_minimal_set.csv',index=False)
    
    transition_df["Minimal Set (same flavour)"] = transition_df["Minimal Set"].apply(rmisc.same_flavour)

    if save_csv == True:
        #transition_df.to_csv('CSV/transition_df_minimal_set_wf.csv',index=False)
        transition_df.to_csv(os.path.join(abcrpv_package_path,"data/main.csv"),index=False)
    print()
    return transition_df

def transition_df():
    try:
        return pd.read_csv(os.path.join(abcrpv_package_path,"data/main.csv"))
    except:
        print("Couldnt't find transition_df in data \nRegenerating",end="")
        return generate_transition_df(save_csv=True)

def get_transitions(mother, daughter):
    return TRANSITION_DF[(TRANSITION_DF["Mother"] == mother) & (TRANSITION_DF["Daughter"]== daughter)]

def generate_transitions_table():
    transitions_table = pd.DataFrame(TRANSITION_DF, columns=['Mother', 'Daughter' , "All Possible Signature",'Minimal Set (same flavour)'])
    try:
        transitions_table['Minimal Set (same flavour)'] = transitions_table['Minimal Set (same flavour)'].map(eval)
        transitions_table['All Possible Signature'] = transitions_table['All Possible Signature'].map(eval)
    except Exception:
        pass
    transitions_table['Minimal Set (same flavour)'] = transitions_table['Minimal Set (same flavour)'].map(frozenset)
    transitions_table['All Possible Signature'] = transitions_table['All Possible Signature'].map(frozenset)
    transitions_table = transitions_table.drop_duplicates(keep='last')
    transitions_table['Minimal Set (same flavour)'] = transitions_table['Minimal Set (same flavour)'].map(list)
    transitions_table['All Possible Signature'] = transitions_table['All Possible Signature'].map(list)
    transitions_table.to_csv(os.path.join(abcrpv_package_path,"data/transitions_table.csv"),index=False)
    return transitions_table

def transitions_table():
    try:
        return pd.read_csv(os.path.join(abcrpv_package_path,"data/transitions_table.csv"))
    except:
        print("Couldnt't find transitions_table in data \nRegenerating...")
        return generate_transitions_table()

def generate_LSP_RPV_decay_table(RPVcoupling): 
# @TODO generate from single to multiple possibilities
    index_name_cat = ['Category', 'LSP', 'decays via', 'Signatures', 'Chain','NV_cascade']
    XSTATE = rdef.STATE_DICT[RPVcoupling.upper()]
    lsp_dec_dat = []
    # j = LSP
    # k = STATE of CAT
    # k[m] = chosen particle.(STATE of CAT)

    for j in rdef.SPARTICLES:
        print(".",end="")
        for k in XSTATE:
                for l in np.arange(0,3,1):
                    rpvdecayed = []
                    rpvdecayed.append(k[0])
                    rpvdecayed.append(k[1])
                    rpvdecayed.append(k[2])
                    rpvdecayed.pop(l)
                    #skip mass degenrate transitions
                    if "W^0"     == j     and "W^+"     == k[l]   :
                        continue
                    if "W^0"     ==  k[l] and "W^+"     == j      :
                        continue
                    if "H^0"     == j     and "H^+"     == k[l]   :
                        continue
                    if "H^0"     ==  k[l] and "H^+"     == j      :
                        continue
                    if "l"       == j     and "nu"      == k[l]   :
                        continue
                    if "l"       == k[l]  and "nu"      == j      :
                        continue
                    if "b_L"     == j     and "t_L"     == k[l]   :
                        continue
                    if "b_L"     == k[l]  and "t_L"     == j      :
                        continue
                    if "nu_tau"  == j     and "tau_L"   == k[l]   :
                        continue
                    if "nu_tau"  == k[l]  and "tau_L"   == j      :
                        continue
                    if j == k[l]:
                        lsp_dec_dat.append([k[-1],j,k[l],rmisc.signature_ordering("".join(list(map(rmisc.sparticles_to_sig, rpvdecayed)))),j+" -- ["+rmisc.sparticles_to_sig(rpvdecayed[0])+","+rmisc.sparticles_to_sig(rpvdecayed[1])+"]"   ,0])
                    else:
                        for i in list(TRANSITIONS_TABLE.loc[(TRANSITIONS_TABLE['Mother'] == j ) & (TRANSITIONS_TABLE['Daughter'] == k[l]) ]["All Possible Signature"].values[0]):
                            chain =(list(TRANSITION_DF.loc[(TRANSITION_DF['Mother'] == j ) & (TRANSITION_DF['Daughter'] == k[l]) & (TRANSITION_DF['Signatures'] == i) ]["Chain"])[0]+" -- ["+rmisc.sparticles_to_sig(rpvdecayed[0])+","+rmisc.sparticles_to_sig(rpvdecayed[1])+"]")
                            nv = list(TRANSITION_DF.loc[(TRANSITION_DF['Mother'] == j ) & (TRANSITION_DF['Daughter'] == k[l]) & (TRANSITION_DF['Signatures'] == i) ]["Nvertex"])[0]
                            lsp_dec_dat.append([k[-1],j,k[l],rmisc.signature_ordering(i.join(list(map(rmisc.sparticles_to_sig, rpvdecayed)))),str(chain),nv])
    lsp_dec_df = pd.DataFrame(lsp_dec_dat,columns=index_name_cat)
    lsp_dec_df["Signatures (ER)"] = lsp_dec_df["Signatures"].apply(rmisc.easy_read)
    lsp_dec_df.to_csv(os.path.join(abcrpv_package_path,"data/"+RPVcoupling+'_1LSP_RPV_DECAY.csv'),index=False)
    print()
    return lsp_dec_df

def one_LSP_RPV_decay_table(RPVcoupling): 
    try:
        return pd.read_csv(os.path.join(abcrpv_package_path,"data/"+RPVcoupling+'_1LSP_RPV_DECAY.csv'))
    except:
        print("Couldnt't find "+RPVcoupling+"_1LSP_RPV_DECAY.csv in data \nRegenerating",end="")
        return generate_LSP_RPV_decay_table(RPVcoupling)

def generate_LSP_sig_cat_table(RPVcoupling):
    if RPVcoupling.upper() == "LLE":
        LSP_dec_RPV_df = LLE_1LSP_RPV_DECAY_TABLE
        rpv_cat = rdef.LLE_CAT
    if RPVcoupling.upper() == "LQD":
        LSP_dec_RPV_df = LQD_1LSP_RPV_DECAY_TABLE
        rpv_cat = rdef.LQD_CAT
    if RPVcoupling.upper() == "UDD":
        LSP_dec_RPV_df = UDD_1LSP_RPV_DECAY_TABLE
        rpv_cat = rdef.UDD_CAT
    onechain_df = pd.DataFrame()
    onechain_df["LSP"] = rdef.SPARTICLES
    onechain_df
    for i in rpv_cat:
        print(".",end="")
        catdat = []
        for j in rdef.SPARTICLES:
            catdat.append(set(np.array(LSP_dec_RPV_df.loc[(LSP_dec_RPV_df['LSP'] == j)& (LSP_dec_RPV_df['Category'] == i)]['Signatures'])))
        onechain_df[i] = catdat
        onechain_df.to_csv(os.path.join(abcrpv_package_path,"data/"+RPVcoupling+'_1LSP_SIG_CAT.csv'),index=False)
    print()
    return onechain_df 

def LSP_sig_cat_table(RPVcoupling): 
    try:
        lsct = pd.read_csv(os.path.join(abcrpv_package_path,"data/"+RPVcoupling+'_1LSP_SIG_CAT.csv'))
        for i in rdef.CAT_DICT[RPVcoupling]:
            lsct[i] = lsct[i].map(eval) 
        return lsct
    except:
        print("Couldnt't find "+RPVcoupling+"_1LSP_SIG_CAT.csv in data \nRegenerating",end="")
        return generate_LSP_sig_cat_table(RPVcoupling)

def generate_2LSP_RPV_decay_table(rpv_coup):
    warnings.filterwarnings("ignore") #ignore panda.append deprecate warning, @TODO update this, do this without append
    lsp_dec_df      = ONE_LSP_RPV_DECAY_DICT[rpv_coup.upper()]
    lsp_onechain_df = ONE_LSP_SIG_CAT_DICT[rpv_coup.upper()]
    rpv_cat         = rdef.CAT_DICT[rpv_coup.upper()]
 
    output_index = ["CAT","LSP A","LSP B","Signature A","Signature B","Chain A","Chain B","Signatures"]
    output_df = pd.DataFrame(columns=output_index)
    for j in rpv_cat:
        print(".",end="")
        for i in rdef.SPARTICLES_DEGENERACY:
            #all possible signatures from one LSP A decay 
            sigseta        = lsp_onechain_df[lsp_onechain_df["LSP"] == i[0]][j].values[0]
            
            #all possible signatures from two of the same LSP A decay 
            pairAA    = list(itertools.product(sigseta, sigseta))
        
            #Get all relevant data
            datAA = [ [j,i[0],i[0],s[0],s[1],
                     list( lsp_dec_df.loc[(lsp_dec_df['LSP'] == i[0] )  & (lsp_dec_df['Category'] == j) & (lsp_dec_df['Signatures'] == s[0])]["Chain"] ),
                     list( lsp_dec_df.loc[(lsp_dec_df['LSP'] == i[0] )  & (lsp_dec_df['Category'] == j) & (lsp_dec_df['Signatures'] == s[1])]["Chain"] ),
                     rmisc.signature_ordering("".join(s))] for s in pairAA]
            entryAA = pd.DataFrame(datAA, columns=output_index)
            output_df = output_df.append(entryAA, ignore_index = True)        
        
            if len(i) > 1:       
                #all possible signatures from one LSP B decay (SU(2) degenerate partner of LSP A)
                sigsetb        = lsp_onechain_df[lsp_onechain_df["LSP"] == i[1]][j].values[0]
                
                #all possible signatures from (LSP A,LSP B) and  (LSP B,LSP B)decay 
                pairAB    = list(itertools.product(sigseta, sigsetb))
                pairBB    = list(itertools.product(sigsetb, sigsetb))
                
                #Get all relevant data
                datAB = [ [j,i[0],i[1],s[0],s[1],
                         list( lsp_dec_df.loc[(lsp_dec_df['LSP'] == i[0] )  & (lsp_dec_df['Category'] == j) & (lsp_dec_df['Signatures'] == s[0])]["Chain"] ),
                         list( lsp_dec_df.loc[(lsp_dec_df['LSP'] == i[1] )  & (lsp_dec_df['Category'] == j) & (lsp_dec_df['Signatures'] == s[1])]["Chain"] ),
                         rmisc.signature_ordering("".join(s))] for s in pairAB]
                datBB = [ [j,i[1],i[1],s[0],s[1],
                         list( lsp_dec_df.loc[(lsp_dec_df['LSP'] == i[1] )  & (lsp_dec_df['Category'] == j) & (lsp_dec_df['Signatures'] == s[0])]["Chain"] ),
                         list( lsp_dec_df.loc[(lsp_dec_df['LSP'] == i[1] )  & (lsp_dec_df['Category'] == j) & (lsp_dec_df['Signatures'] == s[1])]["Chain"] ),
                         rmisc.signature_ordering("".join(s))] for s in pairBB]
                
                entryAB = pd.DataFrame(datAB, columns=output_index)
                entryBB = pd.DataFrame(datBB, columns=output_index)
                output_df = output_df.append(entryAB, ignore_index = True)        
                output_df = output_df.append(entryBB, ignore_index = True)        
    print()
    output_df["Signature A (ER)"] = output_df["Signature A"].apply(rmisc.easy_read)
    output_df["Signature B (ER)"] = output_df["Signature B"].apply(rmisc.easy_read)
    output_df["Signatures (ER)"] = output_df["Signatures"].apply(rmisc.easy_read)
    #output_df.to_csv("CSV/"+rpv_coup.upper()+"_2LSP_table.csv",index=False)
    output_df.to_csv(os.path.join(abcrpv_package_path,"data/"+rpv_coup+'_2LSP_RPV_DECAY.csv'),index=False)

    warnings.filterwarnings("default")
    return output_df

def two_LSP_RPV_decay_table(RPVcoupling): 
    try:
        return pd.read_csv(os.path.join(abcrpv_package_path,"data/"+RPVcoupling+'_2LSP_RPV_DECAY.csv'))
    except:
        print("Couldnt't find "+RPVcoupling+"_2LSP_RPV_DECAY.csv in data \nRegenerating",end="")
        return generate_2LSP_RPV_decay_table(RPVcoupling)

def generate_2LSP_sig_complete(rpv_coup):
    rpv_coup        = rpv_coup.upper().replace(" ","")
    lsp_onechain_df = ONE_LSP_SIG_CAT_DICT[rpv_coup]
    rpv_cat         = rdef.CAT_DICT[rpv_coup]
    
    to_pd_dat = []
    to_pd_index = ["CAT","LSP","AA","AB","BB","ALL",
                        "AA (Common)","AB (Common)","BB (Common)","ALL (Common)",
                        "AA (Same Flavour)","AB (Same Flavour)","BB (Same Flavour)","ALL (Same Flavour)",
                        ]
    for j in rpv_cat:
        print(j)
        for i in rdef.SPARTICLES_DEGENERACY:
            print(".",end="")
            if len(i) > 1:
                #Signatures of 1 LSP decay
                sigseta                  = lsp_onechain_df[lsp_onechain_df["LSP"] == i[0]][j].values[0]
                sigsetb                  = lsp_onechain_df[lsp_onechain_df["LSP"] == i[1]][j].values[0]
                
                #Signatures of 2 LSP decay
                tempsigsetaa             = set([rmisc.signature_ordering("".join(x)) for x in list(itertools.product(sigseta, sigseta))])
                tempsigsetab             = set([rmisc.signature_ordering("".join(x)) for x in list(itertools.product(sigseta, sigsetb))])
                tempsigsetbb             = set([rmisc.signature_ordering("".join(x)) for x in list(itertools.product(sigsetb, sigsetb))])
                tempsigsetall            = tempsigsetbb.union(tempsigsetaa).union(tempsigsetab)

                #Get the minimal set from the 2LSP decay signatures
                ms_tempsigsetaa          = set(minimal_set_greedy(tempsigsetaa))
                ms_tempsigsetab          = set(minimal_set_greedy(tempsigsetab))
                ms_tempsigsetbb          = set(minimal_set_greedy(tempsigsetbb))
                ms_tempsigsetall         = set(minimal_set_greedy(tempsigsetall))

                #Group same flavours the minimal set from the 2LSP decay signatures if exist
                sf_simp_ms_tempsigsetaa  = rmisc.same_flavour(ms_tempsigsetaa) 
                sf_simp_ms_tempsigsetab  = rmisc.same_flavour(ms_tempsigsetab) 
                sf_simp_ms_tempsigsetbb  = rmisc.same_flavour(ms_tempsigsetbb) 
                sf_simp_ms_tempsigsetall = rmisc.same_flavour(ms_tempsigsetall)     


                to_pd_dat.append([j,i,
                                  ms_tempsigsetaa,ms_tempsigsetab,ms_tempsigsetbb,ms_tempsigsetall,
                                  get_common(ms_tempsigsetaa),get_common(ms_tempsigsetab),get_common(ms_tempsigsetbb),get_common(ms_tempsigsetall),
                                  sf_simp_ms_tempsigsetaa,sf_simp_ms_tempsigsetab,sf_simp_ms_tempsigsetbb,sf_simp_ms_tempsigsetall])
                
            else:
                sigset        = lsp_onechain_df[lsp_onechain_df["LSP"] == i[0]][j].values[0]
                tempsigset    = set([rmisc.signature_ordering("".join(x)) for x in list(itertools.product(sigset, sigset))])
                ms_tempsigset = set(minimal_set_greedy(tempsigset))
                sf_simp_ms_tempsigset   = rmisc.same_flavour(ms_tempsigset) 
                to_pd_dat.append([j,i,
                                  ms_tempsigset,"-","-",ms_tempsigset,
                                  get_common(ms_tempsigset),"-","-",get_common(ms_tempsigset),
                                  sf_simp_ms_tempsigset,"-","-",sf_simp_ms_tempsigset])
        print()
    out_2chain_df = pd.DataFrame(to_pd_dat,columns=to_pd_index)
    out_2chain_df.to_csv(os.path.join(abcrpv_package_path,"data/"+rpv_coup.upper()+'_2LSP_sig_complete.csv'),index=False)
    #out_2chain_df.to_csv("CSV/"+rpv_coup.upper()+"_2LSP_sig_complete.csv",index=False)
    return out_2chain_df

def two_LSP_sig_cat_table(RPVcoupling): 
    try:
        lsct = pd.read_csv(os.path.join(abcrpv_package_path,"data/"+RPVcoupling+'_2LSP_sig_complete.csv'))
        #for i in rdef.CAT_DICT[RPVcoupling]:
        #    lsct[i] = lsct[i].map(eval) 
        #return lsct
    except:
        print("Couldnt't find "+RPVcoupling+"_2LSP_sig_complete.csv in data \nRegenerating",end="")
        return generate_2LSP_sig_complete(RPVcoupling)

def minimal_set_greedy(tempsigset,doublecheck=False):
    #start_time = time.time()
    subset        = (get_subsets(tempsigset))
    #print("subset",time.time()-start_time)
    subset_sets   = [x[1] for x in subset]
    #print("subset_sets",time.time()-start_time)
    minimal_set   = rmisc.set_cover_greedy(tempsigset,subset_sets)
    #print("minimal_set",time.time()-start_time)
    
    minimal_set_simplified= []
    for ss in subset:
        for ms in minimal_set:
            if ss[1] == ms:
                minimal_set_simplified.append(ss[0])
                break
    if doublecheck == True:
        leneach = 0
        for k in set(list(map(rmisc.to_L,list(map(rmisc.to_J,tempsigset))))):
            xsig = set([s for s in tempsigset if ((rmisc.to_L((rmisc.to_J(s))))) == k])
            subset = (get_subsets(xsig))
            subset_sets = [x[1] for x in subset]
            minimal_set = rmisc.set_cover_greedy(set(xsig),subset_sets)
            kminimal_set_simplified= []
            for ss in subset:
                for ms in minimal_set:
                    if ss[1] == ms:
                        kminimal_set_simplified.append(ss[0])
                        break
            leneach = leneach + len(kminimal_set_simplified)
            if not set(kminimal_set_simplified).issubset(set(minimal_set_simplified)):
                print("NO")
        if leneach != len(minimal_set_simplified):
            print("NO")

    return minimal_set_simplified

def get_common(insig):
    #input should be list of signature
    #get unique objects in signature
    #print(insig)
    insig = list(itertools.chain.from_iterable(list(map(set_elements_simplified,insig))))
    #print(insig)
    uniquesig = set(list(itertools.chain(*insig))) 
    #print(uniquesig)    
    dat = [] 
    for i in uniquesig:    
        dat.append([i,[x.count(i) for x in insig]])
    return "".join([x[0]*min(x[1]) for x in dat if all(np.array(x[1]) > 0)])

def minimal_sets_simplified_signatures_greedy(insig):
    if len(insig) <= 1:
        return insig, len(insig)
    #Get minimal sets:
    # 1. get all possible subsets that can provide the signatures (get_subsets)
    # 2. from all possible subsets, find minimal sets required to get all siagntures (setCover))
    #input should be list of signautre
    #print(insig)
    outtemp = get_subsets(insig)
    #print(outtemp,list(np.transpose(outtemp)[1]))
    #print((list(np.transpose(outtemp)[1])[0],insig))
    tempcover = rmisc.set_cover_greedy(insig,list(np.transpose(outtemp)[1])[0])
    #print(tempcover)
    output = []
    noutput = []
    for j in tempcover:
        noutput.append(len(j))
        output.append([x[0] for x in outtemp if x[1] == j][0])
    return output,noutput


TRANSITION_DF = transition_df()
TRANSITIONS_TABLE = transitions_table()
try:
    TRANSITIONS_TABLE['Minimal Set (same flavour)'] = TRANSITIONS_TABLE['Minimal Set (same flavour)'].map(eval)
    TRANSITIONS_TABLE['All Possible Signature'] = TRANSITIONS_TABLE['All Possible Signature'].map(eval)
except Exception:
    pass

LLE_1LSP_RPV_DECAY_TABLE = one_LSP_RPV_decay_table("LLE")
LQD_1LSP_RPV_DECAY_TABLE = one_LSP_RPV_decay_table("LQD")
UDD_1LSP_RPV_DECAY_TABLE = one_LSP_RPV_decay_table("UDD")

LLE_1LSP_SIG_CAT_TABLE = LSP_sig_cat_table("LLE")
LQD_1LSP_SIG_CAT_TABLE = LSP_sig_cat_table("LQD")
UDD_1LSP_SIG_CAT_TABLE = LSP_sig_cat_table("UDD")

ONE_LSP_RPV_DECAY_DICT = { "LLE":LLE_1LSP_RPV_DECAY_TABLE,
                           "LQD":LQD_1LSP_RPV_DECAY_TABLE,
                           "UDD":UDD_1LSP_RPV_DECAY_TABLE,}
ONE_LSP_SIG_CAT_DICT = { "LLE":LLE_1LSP_SIG_CAT_TABLE,
                         "LQD":LQD_1LSP_SIG_CAT_TABLE,
                         "UDD":UDD_1LSP_SIG_CAT_TABLE,}

LLE_2LSP_RPV_DECAY_TABLE = two_LSP_RPV_decay_table("LLE")
LQD_2LSP_RPV_DECAY_TABLE = two_LSP_RPV_decay_table("LQD")
UDD_2LSP_RPV_DECAY_TABLE = two_LSP_RPV_decay_table("UDD")


TWO_LSP_RPV_DECAY_DICT = { "LLE":LLE_2LSP_RPV_DECAY_TABLE,
                           "LQD":LQD_2LSP_RPV_DECAY_TABLE,
                           "UDD":UDD_2LSP_RPV_DECAY_TABLE,}


def find_one_lsp_from_signature(signature,rpv_coup="ALL",category="ALL",filename="",save_results=True,verbose=True):
    assert type(signature) == str, "input signature needs to be str, only taking one signature at a time"
    assert rmisc.check_signature_format(signature), "Check signature format, allowed syntax:{x}".format(x=rdef.FINAL_STATE)
    rpv_coup = rpv_coup.upper()
    category = category.upper()
    if category != "ALL":
        if category.count(" ") < 2:
            if category.count("E") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying LLE. Try:\n{x}".format(x=rdef.CAT_DICT["LLE"]))
            if category.count("Q") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying LQD. Try:\n{x}".format(x=rdef.CAT_DICT["LQD"]))
            if category.count("U") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying UDD. Try:\n{x}".format(x=rdef.CAT_DICT["UDD"]))
        if category not in rdef.CAT_DICT[rpv_coup]:
            raise NameError("\""+category+"\" not a category in "+rpv_coup)
    
    if filename == "":
        filename = "one_lsp_that_decays_to_"+signature+"_"+rpv_coup+"_"+category.replace(" ","")+".csv"

    if save_results ==True:
        if ".csv" not in filename:
            filename=filename+".csv"
        results_path = os.path.join(abcrpv_package_path,"results/"+filename)


    rpv_coup = rpv_coup.replace(" ","")
    if category != "ALL" and rpv_coup == "ALL":
        for i,j in rdef.CAT_DICT.items():
            if category in j:        
                rpv_coup = i

    warnings.filterwarnings("ignore") #ignore panda.append deprecategorye warning, @TODO update this, do this without append
    signature = rmisc.signature_ordering(signature)
    output_pd = pd.DataFrame()
    if "L" in signature or "J"in signature or "3"in signature  :
        signature = (set_elements_simplified(signature))
        if verbose == True:
            print("More than one signature in input :\n",filename)
            print("Will be looking up               :\n",signature)
    else:
        signature = [signature]

    if rpv_coup == "ALL":
        for _,LSP_dec_df in ONE_LSP_RPV_DECAY_DICT.items():
            for k in signature:
                output_pd = output_pd.append(LSP_dec_df.loc[(LSP_dec_df['Signatures'] == k)])

    elif rpv_coup != "ALL":
        LSP_dec_df = ONE_LSP_RPV_DECAY_DICT[rpv_coup.upper()]
        if category == "ALL" :
            for k in signature:
                output_pd = output_pd.append(LSP_dec_df.loc[(LSP_dec_df['Signatures'] == k)])

        elif  category != "ALL" :
            for k in signature:
                output_pd = output_pd.append(LSP_dec_df.loc[(LSP_dec_df['Signatures'] == k) & (LSP_dec_df['Category'] == category)])


    if save_results==True:
        if verbose == True:
            print("Results saved in ",results_path)
        output_pd.to_csv(results_path,index=False)     
    return  output_pd  

def find_signatures_from_one_lsp(lsp,rpv_coup="ALL",category="ALL",filename="",save_results=True,verbose=True):
    assert type(lsp) == str        , "input lsp needs to be str, only taking one lsp at a time"
    assert lsp in rdef.SPARTICLES  , "Check LSP format, allowed syntax:{x}".format(x=rdef.SPARTICLES)
    #lsp = rmisc.check_format(lsp)
    rpv_coup = rpv_coup.upper()
    rpv_coup = rpv_coup.replace(" ","")
    category = category.upper()
    if category != "ALL":
        if category.count(" ") < 2:
            if category.count("E") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying LLE. Try:\n{x}".format(x=rdef.CAT_DICT["LLE"]))
            if category.count("Q") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying LQD. Try:\n{x}".format(x=rdef.CAT_DICT["LQD"]))
            if category.count("U") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying UDD. Try:\n{x}".format(x=rdef.CAT_DICT["UDD"]))
        if category not in rdef.CAT_DICT[rpv_coup]:
            raise NameError("\""+category+"\" not a category in "+rpv_coup)

    if filename == "":
        filename = "signatures_from_one_"+lsp.replace("^","")+"_decay_"+rpv_coup+"_"+category.replace(" ","")+".csv"

    if save_results ==True:
        if ".csv" not in filename:
            filename=filename+".csv"
        results_path = os.path.join(abcrpv_package_path,"results/"+filename)

    if category != "ALL" and rpv_coup == "ALL":
        for i,j in rdef.CAT_DICT.items():
            if category in j:        
                rpv_coup = i

    output_pd = pd.DataFrame()
    warnings.filterwarnings("ignore") #ignore panda.append deprecategorye warning, @TODO update this, do this without append
    
    if rpv_coup == "ALL":
        for _,one_lsp_table in ONE_LSP_RPV_DECAY_DICT.items():
        #for one_lsp_table in [LLE_1lsp_table_df,LQD_1lsp_table_df,UDD_1lsp_table_df]:
            output_pd = output_pd.append(one_lsp_table.loc[(one_lsp_table['LSP'] == lsp)])
      
    elif rpv_coup != "ALL":
        one_lsp_table = ONE_LSP_RPV_DECAY_DICT[rpv_coup.upper()]
        if category == "ALL" :
            output_pd = output_pd.append(one_lsp_table.loc[(one_lsp_table['LSP'] == lsp) ])
        elif  category != "ALL" :
            output_pd = output_pd.append(one_lsp_table.loc[(one_lsp_table['LSP'] == lsp) & (one_lsp_table['Category'] == category)])

    if save_results==True:
        if verbose == True:
            print("Results saved in ",results_path)
        output_pd.to_csv(results_path,index=False)     
    return  output_pd              

def find_two_lsp_from_signature(signature,rpv_coup="ALL",category="ALL",filename="",save_results=True,verbose=True):
    assert type(signature) == str, "input signature needs to be str, only taking one signature at a time"
    rpv_coup = rpv_coup.upper()
    category = category.upper()
    rpv_coup = rpv_coup.replace(" ","")
    if filename == "":
        filename = "two_lsp_from_"+signature+"_"+rpv_coup+"_"+category.replace(" ","")+".csv"


    if category != "ALL":
        if category.count(" ") < 2:
            if category.count("E") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying LLE. Try:\n{x}".format(x=rdef.CAT_DICT["LLE"]))
            if category.count("Q") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying LQD. Try:\n{x}".format(x=rdef.CAT_DICT["LQD"]))
            if category.count("U") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying UDD. Try:\n{x}".format(x=rdef.CAT_DICT["UDD"]))
        if category not in rdef.CAT_DICT[rpv_coup]:
            raise NameError("\""+category+"\" not a category in "+rpv_coup)

    if category != "ALL" and rpv_coup == "ALL":
        for i,j in rdef.CAT_DICT.items():
            if category in j:        
                rpv_coup = i

    warnings.filterwarnings("ignore") #ignore panda.append deprecate warning, @TODO update this, do this without append
    signature = rmisc.signature_ordering(signature)
    output_pd = pd.DataFrame()
    if "L" in signature or "J"in signature or "3"in signature  :
        signature = (set_elements_simplified(signature))
        if verbose == True:
            print("More than one signature in input :\n",filename)
            print("Will be looking up               :\n",signature)
    else:
        signature = [signature]

    if rpv_coup == "ALL":
        for _, two_LSP_table in TWO_LSP_RPV_DECAY_DICT.items():
        #for two_LSP_table in [LLE_2LSP_table,LQD_2LSP_table,UDD_2LSP_table]:
            for k in signature:
                output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['Signatures'] == k)])
                #display(two_LSP_table.loc[(two_LSP_table['Signatures'] == k)])

    elif rpv_coup != "ALL":
        two_LSP_table = TWO_LSP_RPV_DECAY_DICT[rpv_coup]   
        if category == "ALL" :
            for k in signature:
                output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['Signatures'] == k)])
                #two_LSP_table.loc[(two_LSP_table['Signatures'] == k)]

        elif  category != "ALL" :
            for k in signature:
                output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['Signatures'] == k) & (two_LSP_table['CAT'] == category)])
    
    warnings.filterwarnings("default") 
    if save_results==True:

        if ".csv" not in filename:
            filename=filename+".csv"

        results_path = os.path.join(abcrpv_package_path,"results/"+filename)

        if verbose == True:
            print("Results saved in ",results_path)
        output_pd.to_csv(results_path,index=False)      
    return output_pd

def find_signatures_from_two_LSP(lspa,lspb,rpv_coup="ALL",category="ALL",filename="",save_results=True,verbose=True):
    assert type(lspa) == str and type(lspb) == str, "input lspa and lspb needs to be str"
    assert lspa in rdef.SPARTICLES  , "Check lspa format, allowed syntax:{x}".format(x=rdef.SPARTICLES)
    assert lspb in rdef.SPARTICLES  , "Check lspb format, allowed syntax:{x}".format(x=rdef.SPARTICLES)
    rpv_coup = rpv_coup.upper()
    rpv_coup = rpv_coup.replace(" ","")
    category = category.upper()
    if category != "ALL":
        if category.count(" ") < 2:
            if category.count("E") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying LLE. Try:\n{x}".format(x=rdef.CAT_DICT["LLE"]))
            if category.count("Q") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying LQD. Try:\n{x}".format(x=rdef.CAT_DICT["LQD"]))
            if category.count("U") > 0:
                raise NameError("Look up rdef.CAT_DICT for categories' syntax, it seems like youre trying UDD. Try:\n{x}".format(x=rdef.CAT_DICT["UDD"]))
        if category not in rdef.CAT_DICT[rpv_coup]:
            raise NameError("\""+category+"\" not a category in "+rpv_coup)
        
    
    if filename == "":
        filename = "signatures_from_2LSP_"+lspa.replace("^","")+"_"+lspb.replace("^","")+"_decay_"+rpv_coup+"_"+category.replace(" ","")+".csv"
    
    if category != "ALL" and rpv_coup == "ALL":
        for i,j in rdef.CAT_DICT.items():
            if category in j:        
                rpv_coup = i

    output_pd = pd.DataFrame()
    warnings.filterwarnings("ignore") #ignore panda.append deprecategorye warning, @TODO update this, do this without append

    if rpv_coup == "ALL":
        for _, two_LSP_table in TWO_LSP_RPV_DECAY_DICT.items():
        #for two_LSP_table in [LLE_2LSP_table,LQD_2LSP_table,UDD_2LSP_table]:
            output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['LSP A'] == lspa) & (two_LSP_table['LSP B'] == lspb)])
            if lspa != lspb:
                output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['LSP A'] == lspb) & (two_LSP_table['LSP B'] == lspa)])
        
    elif rpv_coup != "ALL":
        two_LSP_table = TWO_LSP_RPV_DECAY_DICT[rpv_coup]
        if category == "ALL" :
            output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['LSP A'] == lspa) & (two_LSP_table['LSP B'] == lspb)])
            if lspa != lspb:
                output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['LSP A'] == lspb) & (two_LSP_table['LSP B'] == lspa)])

        elif  category != "ALL" :
            output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['LSP A'] == lspa) & (two_LSP_table['LSP B'] == lspb)& (two_LSP_table['CAT'] == category)])
            if lspa != lspb:
                output_pd = output_pd.append(two_LSP_table.loc[(two_LSP_table['LSP A'] == lspb) & (two_LSP_table['LSP B'] == lspa)& (two_LSP_table['CAT'] == category)])
    
    if save_results==True:
        if ".csv" not in filename:
            filename=filename+".csv"
        results_path = os.path.join(abcrpv_package_path,"results/"+filename)
        if verbose == True:
            print("Results saved in ",results_path)
        output_pd.to_csv(results_path,index=False)     
    warnings.filterwarnings("default") 

    #output_pd.to_csv("RESULTS/signatures_from_"+LSPNAME+".csv")     
    return output_pd

def sanity_checks():
    print("Checking transitions",end="")
    for m in rdef.SPARTICLES:
        print(".",end="")
        for d in rdef.SPARTICLES:
                if m == d:
                     continue
                if (m == "W^+"    and d== "W^0" ) or    (d == "W^+"    and m== "W^0" ):
                      continue
                if (m == "H^+"    and d== "H^0" ) or    (d == "H^+"    and m== "H^0" ):
                      continue
                if (m == "nu"     and d== "l" ) or      (d == "nu"     and m== "l" ):
                      continue
                if (m == "t_L"    and d== "b_L" ) or    (d == "t_L"    and m== "b_L" ):
                      continue
                if (m == "nu_tau" and d== "tau_L" ) or  (d == "nu_tau" and m== "tau_L" ):
                      continue
                if [transition_sig(m,d,i) for i in ("notsup","sup","strsup")].count(None) != 2:
                      print(m,d,[transition_sig(m,d,i) for i in ("notsup","sup","strsup")].count(None))
                      #for i in ("notsup","sup","strsup"):
                      #  print(m,d,i,transition_sig(m,d,i))
    print()                    

    print("Checking functions")
    print("\tChecking find_one_lsp_from_signature",end="")
    for i in rdef.SPARTICLES:
        print(".",end="")
        for j,k in rdef.CAT_DICT.items():
            for l in k:
                if (len(find_signatures_from_one_lsp(i,j,l,save_results=False))) < 1:
                    raise ValueError(i,j,l)
    print()                    
    
    print("\tChecking find_two_lsp_from_signature",end="")
    for i in rdef.SPARTICLES_DEGENERACY:
        print(".",end="")
        if len(i) == 1:
            lspa= i[0]
            lspb= i[0]
        else:
            lspa= i[0]
            lspb= i[1]
        for j,k in rdef.CAT_DICT.items():
                for l in k:
                    if (len(find_signatures_from_two_LSP(lspa,lspb,j,l,save_results=False))) < 1:
                        raise ValueError(lspa,lspb,j,l)
    print()                    

#from . import df_checks
