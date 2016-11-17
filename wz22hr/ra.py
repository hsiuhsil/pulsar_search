
space = 0.15
#ra_scopes = np.arange(np.amin(this_file['POS'][:-1,0]), np.amax(this_file['POS'][:,0]), space)
#dec_scopes = np.arange(np.amin(this_file['POS'][:-1,1]), np.amax(this_file['POS'][:,1]), space)
ra_scopes = np.arange(0., 360., space)
dec_scopes = np.arange(-90., 90., space)
targets = []
for ii in xrange(len(ra_scopes)):
    for jj in xrange(len(dec_scopes)):
        filename = 'RA_'+str(ra_scopes[ii])+'_DEC_'+str(dec_scopes[jj])
        targets.append([filename, ra_scopes[ii], dec_scopes[jj]])

