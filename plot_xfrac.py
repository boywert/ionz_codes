import numpy

def read_xfrac(file):
    f = open(file,"rb")
    #dummy = numpy.fromfile(f,numpy.int32,1)[0]
    ngrid = numpy.fromfile(f,numpy.int32,3)
    #dummy = numpy.fromfile(f,numpy.int32,1)[0]
    #dummy = numpy.fromfile(f,numpy.int32,1)[0]
    data = numpy.fromfile(f,numpy.float32,ngrid[0]**3)
    #dummy = numpy.fromfile(f,numpy.int32,1)[0]
    return data


zlist = open("/mnt/lustre/scratch/cs390/47Mpc/snap_z3.txt").readlines()
for z in zlist:
    z = z.stripe()
    data = read_xfrac("/mnt/lustre/scratch/cs390/codes/ionz_codes/usexfrac/40000.00/xfrac3d_"+z+".bin")
    sum = numpy.sum(data)/306**3
    print z,sum
