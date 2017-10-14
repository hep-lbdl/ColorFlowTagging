from ROOT import *

myfile = TFile("../Singlet_output/singlet_dump-chnk0.root","r")
#myfile = TFile("../Octet_output/octet_dump-chnk0.root","r")
t = myfile.Get("images")

writefile = open("outfile.txt","w")

print t.GetEntries()

for i in range(t.GetEntries()):
    if (i%1000==0):
        print i,t.GetEntries()
        pass
    t.GetEntry(i)
    writefile.write(str(t.pull1)+" "+str(t.pull2)+" -1. -1. -1. -1. "+str(t.jet_m)+" -1. -1. ")
    for j in range(625):
        writefile.write(str(t.image[j])+" ")
        pass
    writefile.write(str(t.jet_delta_R)+" ")
    writefile.write("\n")
