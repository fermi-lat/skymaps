"""
#### defaults for setting up psf parameters #####
Author: M. Roth
$Header$
"""
import pyfits as pf
import numpy as N
import pylab as p
import scipy.optimize

mf = 1.5
mb = 1.2

"""---------------------------Helper functions-------------------------"""
#psf scale function from handoff response
def scale(e,s,p):
	pa=[0,0]

	#front
	if p==0:
		pa[0] = 0.058
		pa[1] = 0.000377
	#back
	if p==1:
		pa[0] = 0.096
		pa[1] = 0.0013
	rs=[]
	for i in N.arange(0,len(s),1):
		rs.append(s[i]*N.sqrt(pa[0]*pa[0]*(e[i]/100)**(-1.6)+pa[1]*pa[1])*180/N.pi)
	return rs

#calculate chi-squared
def chi(a,b,c,d,sarr,Earr):
	logl=0
	for i in N.arange(0,len(sarr)):
		logl+=(sarr[i]-func(a,b,c,d,Earr[i]))**2/(sarr[i]*0.07)**2
	return logl

#fit function
def func(a,b,c,d,le):
	ssq = a**2+c*(le/100)**(b)+d*(le/100)**(2*b)
	if ssq<0 or c<0:
		return 1e80
	return N.sqrt(ssq)

"""----------------------------------------------------------------------------------------"""


"""-------------------------------FITS table lookup----------------------------------------"""
path=r'F:/glast/caldb/v0r7p1/CALDB/data/glast/lat/bcf/psf/'

#names of front and back CALDB files
ffile=path+'psf_P6_v1_diff_front.fits'
bfile=path+'psf_P6_v1_diff_back.fits'

#check for environment variable
import os
if 'CALDB' in os.environ and ffile=='' and bfile=='':
	ffile = os.path.join(os.environ['CALDB'],'v0r7p1','CALDB','data', 'glast', 'lat', 'bsf', 'psf','psf_P6_v1_diff_front.fits')
	bfile = os.path.join(os.environ['CALDB'],'v0r7p1','CALDB','data', 'glast', 'lat', 'bsf', 'bsf','psf_P6_v1_diff_back.fits')

#open fits files and point to tables
frontfile = pf.open(ffile,mode='update')
backfile = pf.open(bfile,mode='update')
fpsftable = frontfile[1].data
bpsftable = backfile[1].data

#x corresponds to energies, y to cos theta
xit = N.arange(0,18,1)
yit = N.arange(0,8,1)
sf=[];sb=[]
gf=[];gb=[]

#go through energies
for i in xit:
	w=[0,0]
	s=[0,0]
	g=[0,0]
	for j in yit:

		#weight contribution for sigma and gamma based on effective area
		#slopes in cos theta are defined at the top of the file as 'mf' and 'mb' for front/back
		costh=(j-7)/10.-0.5
		wf=mf*costh+1
		if wf<0:
			wf=0
		wb=mb*costh+1
		if wb<0:
			wb=0
		w[0]=w[0]+wf
		w[1]=w[1]+wb
		s[0]=s[0]+wf*fpsftable.field('SIGMA')[0][18*j+i]
		s[1]=s[1]+wb*bpsftable.field('SIGMA')[0][18*j+i]

		#weight gamma by core photon fraction
		ncore = fpsftable.field('NCORE')[0][18*j+i]
		if ncore>1:
			ncore=1
		if ncore<0:
			ncore=0
		g[0]=g[0]+wf*ncore*fpsftable.field('GCORE')[0][18*j+i]
		g[0]=g[0]+wf*(1-ncore)*fpsftable.field('GTAIL')[0][18*j+i]
		ncore = bpsftable.field('NCORE')[0][18*j+i]
		if ncore>1:
			ncore=1
		if ncore<0:
			ncore=0
		g[1]=g[1]+wb*ncore*bpsftable.field('GCORE')[0][18*j+i]
		g[1]=g[1]+wb*(1-ncore)*bpsftable.field('GTAIL')[0][18*j+i]
	sf.append(s[0]/w[0])
	sb.append(s[1]/w[1])
	gf.append(g[0]/w[0])
	gb.append(g[1]/w[1])

#setup energy bins as 10/decade starting with 17 MeV (same as CALDB)
energy = 10**(1.25+0.25*xit)

#scale sigma values using PSF scale function from handoff-response
sf = scale(energy,sf,0)
sb = scale(energy,sb,1)

#fit sigma for IParams class
fparams = N.array(scipy.optimize.fmin_powell(lambda x: chi(x[0],x[1],x[2],x[3],sf,energy),(0.01,-1.43,1.46,2.43), ftol = 0.0001,full_output=1, disp=0)[0])
bparams = N.array(scipy.optimize.fmin_powell(lambda x: chi(x[0],x[1],x[2],x[3],sb,energy),(0.01,-1.43,1.46,2.43), ftol = 0.0001,full_output=1, disp=0)[0])

#set gammas
fgam = gf
bgam = gb