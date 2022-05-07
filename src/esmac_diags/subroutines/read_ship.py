#%% read in marmet data
# filename = '../../data/MAGIC/obs/ship/raynolds-marmet/marmet03.txt'
def read_marmet(filename):
    """
    read in marmet data
    """

    f = open(filename, 'r')
    # read in data:
    h = f.readline()
    h = f.readline()
    h = h.strip()
    varlist = h.split()
    data = []
    # if 'data2' in locals():
    #     del(data2)
    for line in f:
        line = line.strip()  # remove \n
        columns = line.split()
        data.append(columns)
        # import numpy as np
        # source  =  []
        # for i in range(0,7):
        #     source.append(columns[i])
        # for i in range(7,len(columns)):
        #     source.append(float(columns[i]))
        # data.append(source)
        # if 'data2' in locals():
        #     data2 = np.column_stack((data2,source))
        # else:
        #     data2 = np.asarray(source)

    f.close()
    # data2[data2< = -999] = np.nan
    return(data, varlist)
