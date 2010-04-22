import os, os.path
import cPickle as pickle
def convert_pickles(picklesDir='../corrections/'):
    """
    NAME:
       convert_pickles
    PURPOSE:
       convert the pickles from containing numpy arrays to containing 
       lists of lists of floats
    INPUT:
       picklesDir - all of the pickles in this directory will be converted
    OUTPUT:
       updated pickles
    HISTORY:
       2010-04-22 - Written - Bovy (NYU)
    """
    for file in os.listdir(picklesDir):
        thisfile= os.path.join(picklesDir,file)
        infile= open(thisfile,'r')
        corr= pickle.load(infile)
        infile.close()
        os.remove(thisfile)
        corrlist= []
        for arr in list(corr):
            corrlist.append([float(a) for a in arr])
        outfile= open(thisfile,'wb')
        pickle.dump(corrlist,outfile)
        outfile.close()

if __name__ == '__main__':
    convert_pickles()
