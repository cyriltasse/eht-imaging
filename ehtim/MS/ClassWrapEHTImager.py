import numpy as np
import ehtim as EHT
import DDFacet.Data.ClassMS
import pyfits
from DDFacet.Other import MyLogger
log = MyLogger.getLogger("ClosureImager")

class ClassWrapEHTImager():
    def __init__(self,**kwargs):
        for key, value in kwargs.items(): setattr(self, key, value)
        self.ReadVisData()
        self.ReadFITSImage()
        self.MakeOBS()

    def ReadVisData(self):
        print>>log,"Reading MS metadata"
        DicoSelectOptions={"FlagAnts":self.FlagAnts.split(",")}
        self.MS = DDFacet.Data.ClassMS.ClassMS(self.MSName, 
                                               Col=self.ColName,
                                               DoReadData=False,
                                               DicoSelectOptions=DicoSelectOptions,
                                               AverageTimeFreq=None)
        print self.MS

        
        DATA = {}
        DATA["iMS"]    = 0
        DATA["iChunk"] = 0

        print>>log,"Reading MS visibilities"
        self.MS.GiveChunk(DATA, 0)

        print>>log,"Selected stations:"
        print>>log,"  %s"%str(self.MS.SelectedStationNames)

        self.TObs=self.MS.DTs+self.MS.dt
        print>>log,"Generating array info (%i stations)"%len(self.MS.SelectedStationNames)
        ArrayTXTName='%s.ArrayPos.txt'%self.MSName
        with open(ArrayTXTName,'w') as f:
            f.write('#NAME X   Y  Z  SEFDR SEFDL FR_PAR_ANGLE FR_ELEV_ANGLE FR_OFFSET[d] DR_RE   DR_IM   DL_RE    DL_IM \n')
            for iAnt,AntName in enumerate(self.MS.SelectedStationNames):
                x,y,z=self.MS.SelectedStationPos[iAnt]
                f.write('%s %f %f %f   1.0001 1.0001 1 -1 0 0 0 0 0\n' % (AntName, x, y, z))

        #self.EHT_Array = EHT.array.load_txt(ArrayTXTName)
        self.EHT_Array = EHT.array.load_txt("LOFAR4.txt")




    def ReadFITSImage(self):
        print>>log,"Read fits image %s"%self.FITSName
        self.EHT_im_prior = EHT.image.load_fits(self.FITSName)

    def MakeOBS(self):
        print>>log,"Generate observation"
        self.EHT_obs = self.EHT_im_prior.observe(self.EHT_Array,10.0,100.0,0.0,12.0,1.0E6, add_th_noise=False)
        #self.TObs=self.MS.DTs+self.MS.dt
        #tint, tadv, tstart, tstop, bw,

        u=self.EHT_obs.data["u"]
        v=self.EHT_obs.data["v"]
        import pylab
        pylab.clf()
        pylab.scatter(u,v)
        pylab.draw()
        pylab.show()

    def main(self):
        map1 = EHT.imager_func(self.EHT_obs,self.EHT_im_prior,self.EHT_im_prior,1.0,d1='amp',d2='cphase',s1='gs',maxit=1000)
        res = self.EHT_obs.res()
        map1blur = map1.blur_gauss((res, res, 0.0),0.7)
        map2 = EHT.imager_func(self.EHT_obs,map1blur,map1blur,1.0,d1='amp',d2='cphase',s1='gs',maxit=300)
        map2blur = map2.blur_gauss((res, res, 0.0),0.5)



#     def Read


# # truth image as sum of Gaussians

# CPTS = np.array([[0.0, 0.0, 10.0, 1.0, 1.0,0.0],\
#                  [-5.0,0.0, 10.5, 0.5, 1.0,0.0],\
#                  [-3.0,0.0,102.5, 2.3, 0.5,90.0],\
#                  [2.0, 0.0, 52.5, 1.0, 1.0,0.0]])

# # prior image - 1 Gaussian in vaguely correct direction

# GPTS = np.array([[-2.0, 0.0, 200.0, 4.0, 0.5,90.0]])
# FOV=100

# ######## utility routine for making gaussians - skip to ### for main parts ############

# def mkgauss (naxes,pos,flux,fwhm,axrat=1.0,angle=0.0,ignore=4.0,dodist=False):
# # note that total flux = peak flux in a pixel * 1.1331*FWHM**2
# # angle is major axis East of North
#     a = np.zeros (naxes[0]*naxes[1]).reshape(naxes[1],naxes[0])
#     fwhm /= 1.66667
#     if axrat==1.0 and angle==0.0:
#         for i in range (naxes[1]):
#             ydist=float(i)-pos[1]
#             for j in range (naxes[0]):
#                 xdist=float(j)-pos[0]
#                 if xdist*xdist+ydist*ydist>ignore*ignore*fwhm*fwhm:
#                     continue
#                 if not dodist:
#                     a[i,j] = flux*np.exp(-(xdist*xdist+ydist*ydist)/ \
#                                     (fwhm*fwhm))/(fwhm*fwhm*np.pi)
#                 else:
#                     a[i,j] = np.hypot(xdist,ydist)
#         return a
#     sinth = np.sin(angle*np.pi/180.0)
#     costh = np.cos(angle*np.pi/180.0)
#     r = np.array([-sinth,costh,-costh,-sinth])
#     rt = np.array([-sinth,-costh,costh,-sinth])
#     sig = np.array([fwhm,0.0,0.0,fwhm*axrat])
#     scr1 = mxmul (sig,r)
#     scr2 = mxmul (rt, scr1)
#     scr1 = mxinv (scr2)
#     for i in range(naxes[1]):
#         ydist=float(i)-pos[1]
#         if abs(ydist)>ignore*fwhm:
#             continue
#         for j in range (naxes[0]):
#             xdist = float(j) - pos[0]
#             if abs(xdist)>ignore*fwhm:
#                 continue
#             ex = scr1[0]*xdist+scr1[1]*ydist
#             ey = scr1[2]*xdist+scr1[3]*ydist
#             if not dodist:
#                 a[i,j] = (flux/axrat)*np.exp(-(ex*ex+ey*ey))/(fwhm*fwhm*np.pi)
#             else:
#                 a[i,j] = np.hypot(ex,ey)/1.6666667
#     return a

# def mxmul(a,b):
#     output=np.zeros(4)
#     output[0]=a[0]*b[0]+a[1]*b[2]
#     output[1]=a[0]*b[1]+a[1]*b[3]
#     output[2]=a[2]*b[0]+a[3]*b[2]
#     output[3]=a[2]*b[1]+a[3]*b[3]
#     return output

# def mxinv(a):
#     det=a[0]*a[3]-a[1]*a[2]
#     output=np.array([a[3],-a[1],-a[2],a[0]])/det
#     return output

# ###########  make image as EHTIM text file from a set of Gaussians #######

# def mkimg(cpts=CPTS,output='simulate.txt',name='LOFAR',ra='00 h 00 m 00.0000 s',\
#           dec='+55 deg 00 m 00.0 s', mjd=48277.000,rf=0.150,fov=100):
#     bigdist = 4.0*np.hypot(cpts[:,0],cpts[:,1]).max()
#     dist = max (bigdist, 4.*cpts[:,3].max())
#     arcperpix = dist/float(fov)
#     a = np.zeros((fov,fov),dtype='float')
#     for c in cpts:
#         x = 0.5*fov+c[0]/arcperpix
#         y = 0.5*fov+c[1]/arcperpix
#         print ('[%d,%d],[%d,%d],%f,%f,%f,%f' % (fov,fov,x,y,c[2],c[3]/arcperpix,c[4],c[5]))
#         a += mkgauss ([fov,fov],[x,y],c[2],c[3]/arcperpix,c[4],c[5])
#     f = open(output,'w')
#     f.write('# SRC: %s\n' % name)
#     f.write('# RA %s\n' % ra)
#     f.write('# DEC %s\n' % dec)
#     f.write('# MJD: %f\n' % mjd)
#     f.write('# RF: %f GHz\n' % rf)
#     f.write('# FOVX: %d pix %f as\n' % (fov,dist))
#     f.write('# FOVY: %d pix %f as\n' % (fov,dist))
#     f.write('# ------------------------------------\n')
#     f.write('# x (as)     y (as)       I (Jy/pixel)  Q (Jy/pixel)  U (Jy/pixel)\n')
#     for i in range(fov):
#         for j in range(fov):
#             f.write('%f %f %f 0.0 0.0\n' % \
#                    (-0.5*dist+i*dist/fov,\
#                     -0.5*dist+j*dist/fov, a[i,j]))
#     f.close()
#     return dist

# ##### make Gaussian prior structure #############

# def mkprior(obs,fov,field,cpts):
#     gaussprior = eh.image.make_square(obs, fov, field/206265.)
#     for cpt in cpts:
#         x, y, fwhm_maj = -cpt[0]/206265., cpt[1]/206265., cpt[3]/206265.
#         fwhm_min = fwhm_maj * cpt[4]
#         theta = np.deg2rad(cpt[5])
#         gaussprior = gaussprior.add_gauss(cpt[2], [fwhm_maj,fwhm_min,theta,x,y])
#     return gaussprior

    
# def test():
#     # make image and uv dataset from CPTS
    
#     field = mkimg (CPTS,fov=FOV)
#     im = eh.image.load_txt('simulate.txt')
#     im.display ()
#     earray = eh.array.load_txt('LOFAR4.txt')

#     obs = im.observe (earray,10.0,100.0,0.0,12.0,1.0E6, add_th_noise=False)
#     stop
#     #obs.plot_bl('ST001','DE601','phase')
#     #obs.plot_cphase('ST001','DE601','DE609')
    
#     # make Gaussian prior image from GPTS
    
#     gaussprior = mkprior (obs,FOV,field,GPTS)
#     gaussprior.display()
    
#     # Two iterations of amp, cphase fitting: blur the output model from 1st and reinsert as prior
#     # need to try some different regularizers here
#     # also need to try different things to minimize in d1 and d2 (and different weights)

#     stop
#     map1 = eh.imager_func(obs,gaussprior,gaussprior,1.0,d1='amp',d2='cphase',s1='gs',maxit=1000)
#     res = obs.res()
#     map1blur = map1.blur_gauss((res, res, 0.0),0.7)
#     map2 = eh.imager_func(obs,map1blur,map1blur,1.0,d1='amp',d2='cphase',s1='gs',maxit=300)
#     map2blur = map2.blur_gauss((res, res, 0.0),0.5)


