import rpv_try.rpv_definitions as rdef
import rpv_try.rpv_misc as rmisc
import rpv_try.abcrpv as rpv

print("Checking transitions",end="")
for m in rdef.SPARTICLES:
    print(".",end="")
    for d in rdef.SPARTICLES:
            if [rpv.transition_sig(m,d,i) for i in ("notsup","sup","strsup")].count(None) < 2:
                  for i in ("notsup","sup","strsup"):
                    print(m,d,i)
print()                    