from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astropy.table import Table,vstack

t1 = ascii.read('NGC1333_Luhman.mrt', format='mrt')
t2 = ascii.read('NGC1333_scholz2012.txt',format='tab')
#we will want to save the extinction later, but the two tables have it in different wavelengths.
#Therefore I create a new column with no values in each table to represent the "other" extincion
t1['AV']=np.nan #add a column with no values to t1
t2['AJ']=np.nan #add a column with no values to t2

print(t1.columns)
print(t2.columns)

#convert the coordinate columns of t1 into format hh:mm:ss dd:mm:ss (similar to t2)
ra1=[str(t1['RAh'][i]) + ':' +  str(t1['RAm'][i]) + ':' +  str(t1['RAs'][i]) for i in range(len(t1))]
dec1=[t1['DE-'][i] + str(t1['DEd'][i]) + ':' +  str(t1['DEm'][i]) + ':' +  str(t1['DEs'][i]) for i in range(len(t1))]

#Convert ra and dec for the two catalogs into format good for matching
c1 = SkyCoord(ra=ra1, dec=dec1, unit=(u.hourangle, u.deg))
c2 = SkyCoord(ra=t2['alpha(J2000)'], dec=t2['delta(J2000)'], unit=(u.hourangle, u.deg))

#We want to save all the common objects + un-matched objects from each of the tables
#match the catalog2 (t2) to catalog1 (t1)
idx, d2d, d3d = c2.match_to_catalog_sky(c1)
max_sep = 1.0 * u.arcsec #matchng radius
idx_sep = d2d < max_sep #all sources that have separation below max_sep arcsec
#in this syntax:
#t2[idx_sep]        -> all matched from the table t2
#t1[idx[idx_sep]]   -> all matched from the table t1
#t2[~idx_sep]       -> all un-matched from the table t2
t1_common=t1[idx[idx_sep]]
t2_only=t2[~idx_sep]
c1_common=c1[idx[idx_sep]] #to keep the coordinates in this format
c2_only=c2[~idx_sep] #to keep the coordinates in this format
t1_common['AV'] = t2[idx_sep]['A_V'] #for the matched elements, update the AV column in t1

#now repeat the same but vice-versa to get t1_only
idx, d2d, d3d = c1.match_to_catalog_sky(c2)
max_sep = 1.0 * u.arcsec #matchng radius
idx_sep = d2d < max_sep #all sources that have separation below max_sep arcsec
t1_only=t1[~idx_sep]
c1_only=c1[~idx_sep]

print('number of objects in the first catalog:', len(t1))
print('number of objects in the second catalog:', len(t2))
print('number of common objects:', np.sum(idx_sep))

#create a new table to save the results
t = Table()
#saving common sources, non-matched from the first, non-matched from the second

ra=np.hstack((c1_common.ra, c1_only.ra, c2_only.ra))
dec=np.hstack((c1_common.dec, c1_only.dec, c2_only.dec))
spt=np.hstack((t1_common['Adopt'], t1_only['Adopt'],t2_only['SpT']))
j=np.hstack((np.array(t1_common['Jmag']), np.array(t1_only['Jmag']),t2_only['J'])) #np.array added because the first arrays have units  and the last one doesn't, so python complained
h=np.hstack((np.array(t1_common['Hmag']), np.array(t1_only['Hmag']))) #the second table doesn't contain H-band magnitude, but it's at the end so it can be left like this
k=np.hstack((np.array(t1_common['Ksmag']), np.array(t1_only['Ksmag']),t2_only['K'])) #np.array added because the first arrays have units  and the last one doesn't, so python complained
av=np.hstack((t1_common['AV'], t1_only['AV'],t2_only['A_V']))
aj=np.hstack((t1_common['Aj'], t1_only['Aj'],t2_only['AJ']))
t['ra'] = ra
t['dec'] = dec
t['spt'] = spt
t['J'] = j
t['H'] = h
t['K'] = k
#also save the extinction. In the first table that's AJ and in the second one Av. Save both
t['AJ']= aj
t['AV']=av

#save the table
t.write('NGC1333.cat',format='ascii',overwrite=True)
