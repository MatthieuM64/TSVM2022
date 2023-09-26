#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import time
import multiprocessing

fontsize=20
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
matplotlib.use('TkAgg')

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

##################
### PARAMETERS ###
##################

#Physical parameters.
v0=0.5
eta=0.24
rho0=0.5
LX=800
LY=100

#Numerical parameters.
init=2
ran=0
DT=25
tmax=50000

#Movie creation.
NCPU=4
multi=True
movie=False

#Read parameters in command line.
for arg in sys.argv[1:]:
	if "-v0=" in arg:
		v0=float(arg[9:])
	elif "-eta=" in arg:
		eta=float(arg[5:])
	elif "-rho0=" in arg:
		rho0=float(arg[6:])
	elif "-LX=" in arg:
		LX=int(arg[4:])
	elif "-LY=" in arg:
		LY=int(arg[4:])
	elif "-init=" in arg:
		init=int(arg[6:])
	elif "-ran=" in arg:
		ran=int(arg[5:])
	elif "-DT=" in arg:
		DT=float(arg[4:])
	elif "-tmax=" in arg:
		tmax=int(arg[6:])
	elif "-NCPU=" in arg:
		NCPU=int(arg[6:])
	elif "-movie"==arg:
		movie=True
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)
		
if NCPU==1:
	multi=False
elif NCPU<1:
	print("Bad value of NCPU: ",NCPU)
	sys.exit(1)
	
rhomax=4
rhomid=0.5*rhomax

#################
### SNAPSHOTS ###
#################

#Create snapshots (frames of the movie)
def Snapshot(i):
	if not os.path.isfile('snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(v0,eta,rho0,LX,LY,init,ran,i)):
		t=i*DT
		
		x=np.linspace(0,LX,LX)
		y=np.linspace(0,LY,LY)
		RHO=np.loadtxt('data_TSVM_dynamics/TSVM_rho_v0=%.8g_eta=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(v0,eta,rho0,LX,LY,init,ran,t))

		fig=plt.figure(figsize=(8,1.5))
		gs=matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],left=0.07,right=1.09, bottom=0.22,top=0.78,hspace=0.1,wspace=0.1)

		ax=plt.subplot(gs[0,0])
		plt.pcolormesh(x,y,RHO,vmin=-rhomax,vmax=rhomax,rasterized=True,cmap=plt.get_cmap('bwr'))
		cb=plt.colorbar(ticks=[-rhomax,0,rhomax],aspect=5)
		cb.ax.set_yticklabels(['$%.4g$'%rhomax,'$0$','$%.4g$'%rhomax])
		cb.solids.set_rasterized(True)
		
		#plt.axis('equal')
		plt.xlim([0,LX])
		plt.ylim([0,LY])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.5*LY,LY])		
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		plt.text(0,1.2*LY,'$t=%d$'%(i*DT),ha='left',va='center')
		plt.text(LX,1.2*LY,'$v_0=%.8g$, $\eta=%.8g$, $\\rho_0=%.8g$'%(v0,eta,rho0),ha='right',va='center')
		
		plt.savefig('snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(v0,eta,rho0,LX,LY,init,ran,i),dpi=480)
		plt.close()
		
		print('-snap=%d/%d -t=%d -rhomin=%.4g -rhomax=%.4g -tcpu=%d'%(i+1,Nsnap,i*DT,RHO.min(),RHO.max(),time.time()-clock))
		del fig,RHO

#Find datafiles and create the corresponding snapshots
os.system('mkdir -p snapshots')

i=0
ARG=[]
while os.path.isfile('data_TSVM_dynamics/TSVM_rho_v0=%.8g_eta=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(v0,eta,rho0,LX,LY,init,ran,i*DT)) and i*DT<=tmax:
	ARG.append(i)
	i+=1

i=0
SNAP=[]
while os.path.isfile('snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(v0,eta,rho0,LX,LY,init,ran,i)) and i*DT<=tmax:
	SNAP.append(i)
	i+=1
	
Nsnap=len(ARG)
print('%d Snapshots (%d already created)'%(Nsnap,len(SNAP)))

if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG)
	pool.close()
	pool.join()
else:
	for i in ARG:
		Snapshot(i)

#############
### MOVIE ###
#############

if movie:	
	os.system('mkdir -p movies')
	#Create the movie
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%%01d.png -c:v h264 -r 25 -s 3840x720 -crf 30 movies/movie_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_LX=%d_LY=%d_init=%d_ran=%d.mp4'%(v0,eta,rho0,LX,LY,init,ran,v0,eta,rho0,LX,LY,init,ran))

print('OK - time=%d sec'%(time.time()-clock))
