import numpy
import os
#parameters
execfile = "./ionz_main"
nion_list = "nion.list"
omegam = 0.27
omegab = 0.044
omegal = 0.73
hubble_h = 0.7
ngrid = 306
boxsize = 47.0
mass_unit =269244.796603

densdir="/scratch/01937/cs390/data/nc306/"
zlistfile="/scratch/01937/cs390/data/snap_z3.txt"
z2listfile="/scratch/01937/cs390/data/snap_z.txt"

def submit_pbs(pbsfile):
    os.system("sbatch "+pbsfile)
    return
def create_pbs(srcdir,outputdir,summaryfile,pbsfile):
    f = open(pbsfile,"w+")
    header = "#!/bin/bash\n"
    header += "#SBATCH -J %s\n" % (pbsfile)
    header += "#SBATCH -o %s.o\%j\n" % (pbsfile)        
    header += """#SBATCH -N 3             
#SBATCH -n 144             
#SBATCH -p normal      
#SBATCH -t 24:00:00        
#SBATCH -A A-asoz"""

    print >> f, header 
    zf = open(zlistfile,"r")
    z3list = zf.readlines()
    zf.close;

    zf = open(z2listfile,"r")
    z2list = zf.readlines()
    zf.close;

    if len(z3list) != len(z2list):
        print "Error: z2 != z2"

    for i in range(len(z3list)):
        z2 = z2list[i].strip()
        z3 = z3list[i].strip()
        if i == 0:
            prev_z = "-1"
        else:
            prev_z = z3list[i-1].strip()
            denfile = densdir+"/"+z3+"n_all.dat"
            srcfile = srcdir+"/"+z2+".dat"
            print >> f, "echo 'z = "+z3+"'"
            print >> f, 'ibrun tacc_affinity',execfile,nion_list,omegam,omegab,omegal,hubble_h,ngrid,boxsize,denfile,srcfile,z3,prev_z,outputdir,summaryfile,mass_unit
    f.close

def main():
    srcdir="/scratch/01937/cs390/data/CSFR/no_reionization/wmap7/"
    outputdir = "/scratch/01937/cs390/data/CSFR/no_reionization/wmap7/SEMNUM/"
    summaryfile = "/scratch/01937/cs390/data/CSFR/no_reionization/wmap7/SEMNUM/sumfile.sum"
    pbsfile = "/scratch/01937/cs390/data/CSFR/no_reionization/wmap7/SEMNUM/lonestar.pbs"
    os.system("mkdir -p "+outputdir)
    create_pbs(srcdir,outputdir,summaryfile,pbsfile)
    #submit_pbs(pbsfile)
    
    srcdir="/scratch/01937/cs390/data/CSFR/okamoto/wmap7/"
    outputdir = "/scratch/01937/cs390/data/CSFR/okamoto/wmap7/SEMNUM/"
    summaryfile = "/scratch/01937/cs390/data/CSFR/okamoto/wmap7/SEMNUM/sumfile.sum"
    pbsfile = "/scratch/01937/cs390/data/CSFR/okamoto/wmap7/SEMNUM/lonestar.pbs"
    os.system("mkdir -p "+outputdir)
    create_pbs(srcdir,outputdir,summaryfile,pbsfile)
    #submit_pbs(pbsfile)

if __name__=="__main__":
    main()
