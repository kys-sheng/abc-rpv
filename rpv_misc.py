"""Misc functions for user interaction."""
try:
    import abcrpv.rpv_definitions as rdef
except:
    import rpv_definitions as rdef
import numpy as np
import itertools

def check_format(instr):
    if "_" in instr:
        x0 = instr.split("_")[0]
        x1 = instr.split("_")[1]
        if x1 == "l" or x1 == "L" or x1 == "left"or x1 == "Left"or x1 == "lh" or x1 ==  "LH":
            x1 = "L" 
        if x1 == "r" or x1 == "R" or x1 == "right" or x1 == "Right" or x1 == "rh" or x1 == "RH":
            x1 = "R"
        if x1.lower() == "tau" or x1 == r"\tau":
            x1 = "tau"
        return x0+"_"+x1
    elif "^" in instr:
        x0 = instr.split("^")[0]
        x1 = instr.split("^")[1]
        if x1 == "+" or x1 == "-" or x1 == "pm" or x1 == r"\pm" or x1.lower() == r"plus" or x1.lower() == r"minus":
            x1 = "pm"
        return x0+"^"+x1
    else:
        return instr

def check_signature_format(instr):
    assert type(instr) == str
    if instr.count(" ") != 0 :
       instr =  instr.replace(" ","")
    if False in [i in rdef.FINAL_STATE for i in list(instr)]:
        return False
    return True

def signature_ordering(instr):
    """
    Make signatures in ordered form based on rdef.sig_order
    @TODO SF bracket treatment

    Inputs:
    ------
    -instr     (str): signature in one-char syntax. Should not have SF brackets (refer same_flavour function) 

    """
    if type(instr) != str:
        return "input not string!"
    
    if instr.count(" ") != 0 :
       instr =  instr.replace(" ","")

    fixedinstr = instr
    nstr = len(instr)
    outstr = ""
    #rdef.sig_order = ["v","J","3","t","b","j","L","T","l","X"]
    
    
    # If () in string, take it out and put it in front
    # Could be done with re, but annoying to learn re syntax (@TODO)
    if "(" in instr and ")" in instr :
        sfdat = ""
        elsedat = ""
        boolget = False
        for i in instr:
            if i == "(":
                boolget = True

            if boolget == True:
                sfdat = sfdat + i
            else:
                elsedat = elsedat + i

            if i == ")":
                boolget = False

        outstr = outstr+sfdat

    else:
        elsedat = instr

    #Reordering
    for i in rdef.SIG_ORDER:
        for j in elsedat:
            if j == i:
                outstr = outstr + i
    
    #String Length should not change by reordering, unless non-signature-syntax char are in string (If so print out)
    if len(outstr) != nstr:
        print("String length changed!",fixedinstr,instr,outstr)
    
    #Two MET is just MET
    while outstr.count("X") > 1 :
        outstr = outstr.replace("XX", "X",100)
        outstr = signature_ordering(outstr)
    
    return outstr


def to_J(inx):
    assert type(inx) == str, "Input must be string!"
    inx = inx.replace("j", "J")
    inx = inx.replace("t", "J")
    inx = inx.replace("b", "J")
    inx = inx.replace("3", "J")
    return inx
    
def to_j3(inx):
    assert type(inx) == str, "Input must be string!"
    inx = inx.replace("t", "3")
    inx = inx.replace("b", "3")
    return inx
    
def to_L(inx):
    assert type(inx) == str, "Input must be string!"
    inx = inx.replace("l","L")
    inx = inx.replace("T","L")
    return inx

def check_merge(insiglist,merge_choice):
    """
    Check if merge (group signatures via simplified object (J,3,L)) should be considered:
    Conditions:
    1: if in total, at least one of each of oblist is in the whole set
    2: if simplified exist in every element 
    @TODO SF bracket treament

    Inputs
    ------
    -insiglist    (array/list of str) : array/list of signatures in one-char syntax
    -merge_choice               (str) : "J"/"3"/"L"
    
    Output:bool
    """
    oblist = []
    if merge_choice == "J":
        oblist = ["j","t","b"]
    if merge_choice == "3":
        oblist = ["t","b"]
    if merge_choice == "L":
        oblist = ["l","T"]
    nobdat = []
    for i in oblist:
        nob = 0
        for j in insiglist:
            nob = nob + j.count(i)
        nobdat.append(nob)

    return all(x > 0 for x in nobdat) or all (x.count(merge_choice) > 0 for x in insiglist)


def easy_read(instr):
    """
    Go from one-char syntax to easyread syntax
    @TODO SF bracket treatment
    
    Inputs:
    ------
    -instr     (str): signature in one-char syntax. Should not have SF brackets (refer same_flavour function)
    
    Example:
    ------
    jjjJ3TLEv -> 1J + 3j_l + 1j_3 + 1L + 1tau + 1v + MET
    
    """
    output = ""
    orderdat = [instr.count("J"),"J"],[instr.count("j"),"j_l"],[instr.count("3"),"j_3"],[instr.count("t"),"t"],[instr.count("b"),"b"],[instr.count("l"),"l"],[instr.count("L"),"L"],[instr.count("T"),"tau"],[instr.count("v"),"v"],[instr.count("X"),"MET"]
    for i in orderdat:
        if i[1] == "MET" and i[0] == 0:
            output = output[:-2]
        elif i[0] == 0:
            continue
        elif i[1] == "MET" and i[0] != 0:
            output = output + " MET"
        elif i[1] != "MET":
            output = output + " " + str(i[0])  + str(i[1]) + " +" 
    return output

def setCover(setList,target=None):
    """
    Bruteforce method for Set Cover problem (Not ideal when signatures goes up to >9 objects)
    https://stackoverflow.com/questions/60114801/find-minimal-set-of-subsets-that-covers-a-given-set
    
    Inputs
    ------
    -setList    : list of subsets
    
    Output:
    ------
    -bestCover  : minimal list of subsets required to cover whole set
    """    
    if not setList: return None
    if target  is None: target  = set.union(*setList)
    bestCover = []
    for i,values in enumerate(setList):
        remaining = target - values
        if remaining == target: continue
        if not remaining: return [values]
        subCover = setCover(setList[i+1:],remaining)
        if not subCover: continue
        if not bestCover or len(subCover)<len(bestCover)-1:
            bestCover = [values] + subCover
    return bestCover


def same_flavour(insig):
    #input should be list of signature
    output = []
    for iii in set(list(map(to_L,set(list(map(to_J,insig)))))):
        #group signature by same number of objects        
        xsig = [s for s in insig if ((to_L((to_J(s))))) == iii]
        #print("xsig",xsig)
        
        # if only one element in group, nothing to do
        if len(xsig) == 1:
            output.append(xsig[0])
            continue

        #get unique objects in signature
        uniquesig = set(list(itertools.chain(*xsig))) 
        #print("uniquesig",uniquesig)
        
        samedat = ""
        #print("samedat",samedat)
        dat = [] 
        #print("dat",dat)
        
        for i in uniquesig:    
            if all(x.count(i) > 0 for x in xsig) : 
                samedat = samedat + i * int(min(list(set([x.count(i) for x in xsig]))))
                #print("samedat",samedat)
                if len(set([x.count(i) for x in xsig])) != 1:
                    leftdat = np.array([x.count(i) for x in xsig]) - np.array(len(xsig)*[int(min(list(set([x.count(i) for x in xsig]))))])
                    #print("leftdat",leftdat)
                    dat.append([x*i for x in leftdat])
            else:
                 dat.append([x.count(i)*i for x in xsig])
        diffdat = [signature_ordering("".join(x)) for x in np.transpose(dat)]
        #print("diffdat",diffdat)
        #samedat is common objects in group 
        #diffdat is different objects in group
        
        if set(diffdat) == set(["jj","tt","bb"]) or set(diffdat) == set(["jj","33"]):
            diffdat = "(JJ)"
            #print("diffdat",diffdat)
            output.append( diffdat+signature_ordering(samedat))
        elif set(diffdat) == set(["tt","tb","bb"]):
            diffdat = "(33)"
            #print("diffdat",diffdat)
            output.append( diffdat+signature_ordering(samedat))
        #elif set(diffdat) == set(["tt","bb"]):
        #    diffdat = r"{33}"
        #    #print("diffdat",diffdat)
        #    output.append( diffdat+signature_ordering(samedat))
        elif set(diffdat) == set(["ll","TT"]):
            diffdat = "(LL)"
            #print("diffdat",diffdat)
            output.append( diffdat+signature_ordering(samedat))
        else:
            for i in diffdat:
                output.append(i+signature_ordering(samedat)) 
                
    return output

def sparticles_to_sig(instr):
    if instr == "q" or instr == "d" or instr == "u" or instr == "d_L" or instr == "u_L" :
        instr = "j"
    if instr == "b_L" or instr == "b" :
        instr = "b"
    if instr == "t_L" or instr == "t"  :
        instr = "t"
    if instr == "l" or instr == "e" :
        instr = "l"
    if instr == "nu" or instr == "nu_tau":
        instr = "X"
    if instr == "tau_L" or instr == "tau":
        instr = "T"
    return instr 

def set_cover_greedy(universe, subsets):
    """Find a family of subsets that covers the universal set"""
    elements = set(e for s in subsets for e in s)
    # Check the subsets cover the universe
    if elements != universe:
        return None
    covered = set()
    cover = []
    # Greedily add the subsets with the most uncovered points
    while covered != elements:
        subset = max(subsets, key=lambda s: len(s - covered))
        cover.append(subset)
        covered |= subset
 
    return cover

def flatten_list(list_of_lists,flat_list=None,level=False):
    if level == False:
         flat_list=[]
    for item in list_of_lists:
            if type(item) == list:
                flatten_list(item, flat_list,True)
            else:
                flat_list.append(item)

    return flat_list