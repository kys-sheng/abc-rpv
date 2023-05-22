import numpy as np
import pandas as pd  
from IPython.display import display 
import subprocess
import itertools
import warnings
import os 
import rpv_definitions as rdef
import rpv_misc as rmisc
from pathlib import Path

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
                      - data/table_notsup.dat
                      - data/table_strsup.dat
                      - data/table_sup.dat
    -opt            : print message if cant find it(default=False)
    """
    trans_pd = pd.read_csv(os.path.join(Path(__file__).parent.absolute(),"data/table_"+sup+".csv"))
    trans_pd[(trans_pd["Mother"] == mother) & (trans_pd["Daughter"] == daughter)]
    transition = trans_pd[(trans_pd["Mother"] == mother) & (trans_pd["Daughter"] == daughter) ]["Transitions"].values
    if len(transition) == 0:
        if opt == True:
            print("\n"+(rmisc.check_format(mother)+" to "+rmisc.check_format(daughter)+" is not "+sup))
        return None
    assert len(transition) == 1, ("Should be one occurance of the transition in table_{sup}.csv".format(sup=sup), transition)
    return transition[0].split()