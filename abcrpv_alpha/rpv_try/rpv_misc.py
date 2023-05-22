"""Misc functions for user interaction."""
import rpv_definitions as rdef

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

def signature_format_check(instr):
    assert type(instr) == str
    if False in [i in rdef.FINAL_STATE for i in list(instr)]:
        return False
    return True

def signature_ordering(instr):
    """
    Make signatures in ordered form based on sig_order
    @TODO SF bracket treatment

    Inputs:
    ------
    -instr     (str): signature in one-char syntax. Should not have SF brackets (refer same_flavour function) 

    """
    fixedinstr = instr
    nstr = len(instr)
    outstr = ""
    sig_order = ["v","J","3","t","b","j","L","T","l","E"]
    if type(instr) != str:
        return "input not string!"
    
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
    for i in sig_order:
        for j in elsedat:
            if j == i:
                outstr = outstr + i
    
    #String Length should not change by reordering, unless non-signature-syntax char are in string (If so print out)
    if len(outstr) != nstr:
        print("String length changed!",fixedinstr,instr,outstr)
    
    #Two MET is just MET
    while outstr.count("E") > 1 :
        outstr = outstr.replace("EE", "E",100)
        outstr = signature_ordering(outstr)
    
    return outstr

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
    orderdat = [instr.count("J"),"J"],[instr.count("j"),"j_l"],[instr.count("3"),"j_3"],[instr.count("t"),"t"],[instr.count("b"),"b"],[instr.count("l"),"l"],[instr.count("L"),"L"],[instr.count("T"),"tau"],[instr.count("v"),"v"],[instr.count("E"),"MET"]
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


