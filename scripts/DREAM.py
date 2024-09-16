# -*- coding: utf-8 -*-
"""
Created on Tue May  9 20:32:52 2017

@author: James
"""

import numpy as np
import scipy as sp

class DREAM:
    def __init__(self, nchains, npars,nburn = 100, npairs = 1,prior = False, k =5):
        self.nc = nchains
        self.npars = npars
        self.CR = 0.5
        self.chains = []
        self.genr = 0
        self.delta = npairs
        self.nb = nburn
        self.burn = True
        self.ct = 0
        self.k = k
        for i in range(nchains):
            self.chains.append(chain(npars,prior = prior))
            
    def par_set(self,log,uniform,mean,width,pmin,pmax):
        for i in range(self.nc):
            self.chains[i].par_set(log,uniform,mean,width,pmin,pmax,rstart = True)
            
    def gen_mod(self):
        self.genr += 1
        if self.genr == 5:
            self.genr = 0
            
    def par_reset(self):
        for i in range(self.nc):
            self.chains[i].pars = self.chains[i].current 
            
    def init_set(self,ninit):
        ini = np.zeros((ninit,self.npars))
        for j in range(ninit):
            for i in range(self.npars):
                if self.chains[0].uniform[i] == True:
                    ini[j,i] = self.chains[0].mean[i] - self.chains[0].width[i] + np.random.rand() * 2. *self.chains[0].width[i]
                else:
                    ini[j,i] = np.random.normal(loc = self.chains[0].mean[i], scale = self.chains[0].width[i])
        return(ini)
        
    def normal_sampling(mu, sigma, nr_param):
        nums = [] 
        count = 0
        while count < nr_param:
            temp = random.gauss(mu,sigma) 
            if temp >= mu - 3 * sigma and temp <= mu + 3*sigma:
                nums.append(temp)
                count += 1
        return(nums)  

    def propgen(self):
        if self.burn == True:
            self.dumloc = np.random.randint(self.ncr) 
            self.CR = (self.dumloc +1) / self.ncr
        for i in range(self.nc):
            #Identify random chains to use
            dum = list(range(self.nc))
            dum.remove(i)
            dx = np.zeros(self.npars)
            ll = self.nc-1
            if self.delta < 0:
                npairs=np.random.randint(1,abs(self.delta)+1)
            else:
                npairs = self.delta
            for k in range(npairs):
                dum2  = np.random.randint(ll)
                ll = ll-1
                dx += self.chains[dum[dum2]].current
                #print(dum2,dx)
                dum.remove(dum[dum2])
                dum2  = np.random.randint(ll)
                ll = ll-1
                dx += -self.chains[dum[dum2]].current
                #print(dum2,dx)
                dum.remove(dum[dum2])
            #print(dx)
            #dx = dx/npairs
            d = self.npars
            #print(d)
            mult = np.ones(self.npars)
            for k in range(self.npars):
                if np.random.rand() > self.CR:
                    dx[k] = 0.
                    d += -1
                    mult[k] = 0.
            if d<1:
                d = 1
            if self.genr < 4:
                gamma = 2.38/np.sqrt(2.*npairs*d)
            else:
                gamma = 1.0
            if gamma > 1.0:
                gamma = 1.0
            dx = dx * gamma # + np.random.normal(scale = self.chains[i].pwidth) * mult
            #print(gamma,dx)
            self.chains[i].proposal = self.chains[i].current + dx
            for k in range(self.npars):
                while self.chains[i].proposal[k] < self.chains[i].pmin[k] or elf.chains[i].proposal[k] > self.chains[i].pmax[k]:
                    if self.chains[i].proposal[k] < self.chains[i].pmin[k]:
                         self.chains[i].proposal[k] = self.chains[i].pmax[k] - (self.chains[i].pmin[k] - self.chains[i].proposal[k])
                    elif self.chains[i].proposal[k] > self.chains[i].pmax[k]:
                         self.chains[i].proposal[k] = self.chains[i].pmin[k] + (self.chains[i].proposal[k]-self.chains[i].pmax[k])
        self.ct+=1
    
    def propgenzs(self,current,lookback = None):
        if self.burn == True:
            self.dumloc = np.random.randint(self.ncr) 
            self.CR = (self.dumloc +1) / self.ncr
        if self.delta < 0:
            npairs=np.random.randint(1,abs(self.delta)+1)            
        else:
            npairs = self.delta  
        dum = np.zeros((self.k,self.npars))
        if self.genr <4:
            gamma = 2.38/np.sqrt(2.*npairs*self.npars)
        else:
            gamma = 1.
        mean = np.mean(self.chains[0].pars,axis=0)
        for i in range(1,self.nc):
            mean += np.mean(self.chains[i].pars,axis=0)
        mean /= self.nc
        var = np.zeros_like(mean)
        for i in range(self.nc):
            for j in range(self.npars):
                var[j] += np.sum((self.chains[i].pars[:,j]-mean[j])**2)
        var = var/(self.nc*self.npars-1)
        b = np.sqrt(var)/1e3 
        b[b==0] = 1e-5
        if gamma>1.:
            gamma = 1.
        b1 = np.zeros_like(b)
        b2 = np.zeros_like(b)
        for i in range(self.npars):
            b1[i] = -b[i] + np.random.rand()*2*b[i]
            b2[i] =  np.random.normal(0,b[i])
            
        dx = np.zeros(self.npars)
        #b1 = 0.00001
        #b2 = 0.00001
        for k in range(self.k):
            for j in range(npairs):
                ii = np.random.randint(self.nc)
                if lookback == None:
                    lb = np.shape(self.chains[ii].pars)[0]
                elif lookback> np.shape(self.chains[ii].pars)[0]:
                    lb = np.shape(self.chains[ii].pars)[0]
                else:
                    lb = np.copy(lookback)
                #print(lb)
                jj = np.random.randint(lb)
                dx+= self.chains[ii].pars[-jj-1,:]
                ii = np.random.randint(self.nc)
                if lookback == None:
                    lb = np.shape(self.chains[ii].pars)[0]
                elif lookback> np.shape(self.chains[ii].pars)[0]:
                    lb = np.shape(self.chains[ii].pars)[0]
                else:
                    lb = np.copy(lookback)
                #print(lb)
                jj = np.random.randint(lb)
                dx-= self.chains[ii].pars[-jj-1,:]
            for j in range(self.npars):
                if np.random.rand() > self.CR:
                    dx[j] = 0.
            dum[k,:] = current + (1+b1)*gamma*dx+b2
            for j in range(self.npars):
                while True:
                    if dum[k,j] < self.chains[0].pmin[j]:
                         dum[k,j] = self.chains[0].pmax[j] - (self.chains[0].pmin[j] - dum[k,j])
                    elif dum[k,j] > self.chains[0].pmax[j]:
                         dum[k,j] = self.chains[0].pmin[j] + (dum[k,j]-self.chains[0].pmax[j])
                    if dum[k,j] >= self.chains[0].pmin[j] and dum[k,j] <= self.chains[0].pmax[j]:
                        break
        
                    
        return(dum)
            
        
        
    def accept_rate(self):
        np = 0
        nc = 0
        for i in range(self.nc):
            np += self.chains[i].gen
            nc += self.chains[i].trans
        self.rate = nc/np
            
    
    def delm_update(self):
        for i in range(self.nc):
            for j in range(self.npars):
                self.delm[self.dumloc] += (self.chains[i].pars[-1,j] - self.chains[i].pars[-2,j]) **2. /self.Var[j]                
                self.crct[self.dumloc] +=1
                
    def Chain_removal(self):
        self.reset = False
        Omega = np.zeros((self.nc))
        n = int(self.ct/2)
        last = np.zeros((self.nc))
        for i in range(self.nc):
            Omega[i] = np.average(self.chains[i].likelihood[n:])
            last[i] = self.chains[i].likelihood[-1]
        best = np.argmax(Omega)
        IQRmin = np.percentile(Omega,25) - 2 * (np.percentile(Omega,75) - np.percentile(Omega,25))        
        for i in range(self.nc):      
            if Omega[i] < IQRmin:
                print(Omega[i],IQRmin)
                self.chains[i].current = np.copy(self.chains[best].current)
                self.chains[i].Lold = np.copy(self.chains[best].Lold)
                self.reset = True
        if self.reset == False and self.ct > self.nb:
            self.burn = False
        elif self.reset == True:
            self.reset = False
            self.ct = -1
            self.genr = 0
            for i in range(self.nc):
                self.chains[i].likelihood = np.array([self.chains[i].Lold])
                self.chains[i].pars = np.array([self.chains[i].current])
                

    def set_CR(self, ncr = 3):
        self.crct=np.zeros(ncr)
        self.ncr = ncr
        self.delm = np.zeros(ncr)
        
    def Rget(self):
        mean = np.zeros((self.nc,self.npars))
        var = np.zeros((self.nc,self.npars))
        for i in range(self.nc):
            mean[i,:], var[i,:], n = self.chains[i].varget() 
        W = np.average(var,axis=0)
        B = np.var(mean,axis=0) * n
        self.Var = (1. -1./n) * W + 1. /n * B
        self.R = np.sqrt(self.Var/W)
        
        
class chain:
    def __init__(self,npars,prior = False):
        self.npars = npars
        self.current = np.zeros(npars)
        self.proposal = np.zeros(npars)
        self.Lold = 0.
        self.Lnew = 0.
        self.pwidth = 0.
        self.prior = prior
        self.likelihood = []
        self.trans = 0
        self.gen = 0
        
    def par_set(self,log,uniform,mean,width,pmin,pmax,rstart = False):
        self.log = log #logical
        self.uniform = uniform # logical
        self.mean = mean # mean or centre of distribution (log for log parameters)
        self.width = width #stabdard deviation or width (;og for log)
        self.rs = rstart #logica - if you want random generte start point
        self.pwidth = []    
        self.pmin=pmin
        self.pmax=pmax
        for i in range(np.size(self.width)):        
            self.pwidth.append(self.width[i]/25.)
        if self.rs == True:
            for i in range(self.npars):
                if uniform[i] == True:
                    self.current[i] = self.mean[i] - self.width[i] + np.random.rand() * 2. *self.width[i]
                else:
                    self.current[i] = np.random.normal(loc = self.mean[i], scale = self.width[i])
        else:
            self.current[i] = np.copy(self.mean[i])
        for i in range(self.npars):
            if self.current[i] < self.pmin[i]:
                self.current[i] = 2. * self.pmin[i] - self.current[i]
            if self.current[i] > self.pmax[i]:
                self.current[i] = 2. * self.pmax[i] - self.current[i]                
        self.pars = np.zeros((1,self.npars))
        self.pars[0,:] = self.current
        
        
        
    def prior(self):
        for i in range(self.npars):
            if self.uniform == True:
                if self.proposal[i] < (self.mean[i] - self.width[i]) or self.proposal[i] > (self.mean[i] + self.width[i]):
                    self.Lnew += -1000
            else:
                self.Lnew += ((self.proposal[i]-self.mean[i])**2.)/(2.*self.width[i]**2.)
                
    def tprob(self):
        self.gen += 1
        if self.Lnew > self.Lold:
            self.current = np.copy(self.proposal)
            self.Lold = np.copy (self.Lnew)
            self.trans+=1
        else:
            if np.random.rand() < np.exp(self.Lnew - self.Lold):
                self.current = np.copy(self.proposal)
                self.Lold = np.copy(self.Lnew)
                self.trans+=1
        self.pars = np.vstack((self.pars,self.current))
        self.likelihood.append(self.Lold)
        
    def varget(self):
        n = int(np.shape(self.pars)[0]/2)
        s2 = np.var(self.pars[n:,:],axis = 0)
        mn = np.average(self.pars[n:,:], axis=0)
        return(mn,s2,n)
            
                
        