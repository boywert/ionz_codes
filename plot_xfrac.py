import numpy

def read_xfrac(file):
    f = open(file,"rb")
    #dummy = numpy.fromfile(f,numpy.int32,1)[0]
    ngrid = numpy.fromfile(f,dtype=numpy.int32,3)
    #dummy = numpy.fromfile(f,numpy.int32,1)[0]
    #dummy = numpy.fromfile(f,numpy.int32,1)[0]
    data = numpy.fromfile(f,dtype=numpy.float32,ngrid*ngrid*ngrid)
    #dummy = numpy.fromfile(f,numpy.int32,1)[0]
    return data

data = read_xfrac("/mnt/lustre/scratch/cs390/codes/ionz_codes/usexfrac/40000.00/xfrac3d_6.354.bin")

sum = numpy.sum(data)/306**3
print sum
