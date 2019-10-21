#!/usr/bin/env python
# coding: utf-8

#Model of expression with an inductor at 3uM
# In[1]:


def pend(z, t, K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, Ym, n, PoPS, Ymin, Ymax, K, n2, a, b, c, Ecr, EEcr, EEEcr, j_m_V, j_m_P, j_m_C, j_m_M, Vo, Po, Co, Mo):
    Vc, Pc, Cc, Mc, R1, Ec, EEc, EEEc = z #Variables to integrate
    dVcdt = -j_m_V*(Vc-Vo)/volc-(K1cat*Ec*Vc)/(Vc + K1m) #Vanilic acid
    dPcdt = -j_m_P*(Pc-Po)/volc+(K1cat*Ec*Vc)/(Vc + K1m) - (K2cat*EEc*Pc)/(Pc + K2m) #PCA
    dCcdt = -j_m_C*(Cc-Co)/volc+(K2cat*EEc*Pc)/(Pc + K2m) - (K1cat*EEEc*Cc)/(Cc + K3m) #Catechol
    dMcdt = -j_m_M*(Mc-Mo)/volc+(K1cat*EEEc*Cc)/(Cc + K3m) #Cis,Cis-Muconic Acid
    dR1dt = (n*PoPS*(Ymin+(Ymax-Ymin)*(((Vc*1000)**n2)/((K)**n2+(Vc*1000)**n2))))/(volc*6.022e20)-Ym*R1 #mRNA
    dEcdt = a*R1-Ec*Ecr #ligM monomer
    dEEcdt = b*R1-EEc*EEcr #aroY monomer
    dEEEcdt = c*R1-EEEc*EEEcr #catA monomer
    dCdt = [dVcdt,dPcdt,dCcdt,dMcdt,dR1dt,dEcdt,dEEcdt,dEEEcdt] #Functions to integrate
    return dCdt


# In[2]:


K1cat = (0.078*9.14) #Kcat LigM
K1m = 0.078 #Km LigM
K2cat = 0.00009 #Kcat AroY
K2m = 0.59 #Km AroY
K3cat = (0.00185*10485) #Kcat CatA
K3m = 0.00185 #Kmm CatA
volo = 1000000 #Medium volume (outside)
volc = 1e-15 #Citoplasm volume
volb = 0 #BMC volume
Ym = 0.0043 #mRNA degradation rate
n = 15 #Number of plasmid copies
PoPS = 0.03 #Transcription rate
Ymin = 2.4e-3 #Min transcription rate (RPU's)
Ymax = 3 #Max transcription rate (RPU's)
K = 26 #Marionette constant (uM)
n2 = 2.3 #Another marionette constant
a = 29.97 #ligM traduction rate
b = 26.30 #aroY traduction rate
c = 41.30 #catA traduction rate
Ecr = 0.000385 #ligM dilution rate
EEcr = 0.000385 #aroY dilution rate
EEEcr = 0.000385 #catA dilution rate
Vo = 3e-3 #vanilic acid outside
Po = 0 #protocatechuic acid outside
Co = 0 #catechol outside
Mo = 0 #cis,cis-muconic acid outside
#flux constants for species trough membrane
j_m_V = 1e-13 # Va
j_m_P = 1e-13 # PCA
j_m_C = 1e-13 # Catechol
j_m_M = 1e-13 #MA


# In[3]:


z0 = [0,0,0,0,0,0,0,0] #Initial values


# In[4]:


import numpy as np
t = np.linspace(0, 0.1, 100001) #Create 100001 values between 0-0.1 in the X axis


# In[5]:


from scipy.integrate import odeint #Integrate all functions for the X axis given values
sol1 = odeint(pend, z0, t, args=(K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, Ym, n, PoPS, Ymin, Ymax, K, n2, a, b, c, Ecr, EEcr, EEEcr, j_m_V, j_m_P, j_m_C, j_m_M, Vo, Po, Co, Mo))


# In[6]:


z0 = list(sol1[100000,:]) #New initial values using final values for past function


# In[7]:


t = np.linspace(0.1, 10, 100001) #Give new values to the X axis for the next step


# In[8]:


from scipy.integrate import odeint
sol2 = odeint(pend, z0, t, args=(K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, Ym, n, PoPS, Ymin, Ymax, K, n2, a, b, c, Ecr, EEcr, EEEcr, j_m_V, j_m_P, j_m_C, j_m_M, Vo, Po, Co, Mo))


# In[9]:


z0 = list(sol2[100000,:])


# In[10]:


t = np.linspace(10, 1000, 100001)


# In[11]:


from scipy.integrate import odeint
sol3 = odeint(pend, z0, t, args=(K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, Ym, n, PoPS, Ymin, Ymax, K, n2, a, b, c, Ecr, EEcr, EEEcr, j_m_V, j_m_P, j_m_C, j_m_M, Vo, Po, Co, Mo))


# In[12]:


z0 = list(sol3[100000,:])


# In[13]:


t = np.linspace(1000, 14400, 100001)


# In[14]:


from scipy.integrate import odeint
sol4 = odeint(pend, z0, t, args=(K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, Ym, n, PoPS, Ymin, Ymax, K, n2, a, b, c, Ecr, EEcr, EEEcr, j_m_V, j_m_P, j_m_C, j_m_M, Vo, Po, Co, Mo))


# In[15]:


solprueba= np.concatenate((sol1,sol2), axis=0) #This concatenate functions are used for mix the four parts in one single variable


# In[16]:


solprueba2= np.concatenate((sol3,sol4), axis=0)


# In[17]:


sol= np.concatenate((solprueba,solprueba2), axis=0) 


# In[18]:


a = np.linspace(0, 0.1, 100001)
b = np.linspace(0.1, 10, 100001)
c = np.linspace(10, 1000, 100001)
d = np.linspace(1000, 14400, 100001)
e = np.concatenate((a, b), axis=None)
f = np.concatenate((c, d), axis=None)
t = np.concatenate((e, f), axis=None) #All values for X axis


# In[19]:


import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
plt.plot(t, sol[:, 0], 'b-', label='V.A.-cell')
plt.plot(t, sol[:, 1], 'g-', label='PCA-cell')
plt.plot(t, sol[:, 2], 'r-', label='C-cell')
plt.plot(t, sol[:, 3], 'm-', label='M.A.-cell')
plt.plot(t, sol[:, 4], 'b--', label='mRNA')
plt.legend(loc='best')
plt.xlabel('t(s)')
plt.ylabel('C(mM)')
plt.grid()
plt.show()


# In[20]:


plt.plot(t, sol[:, 1], 'g-', label='PCA-cell')
plt.plot(t, sol[:, 2], 'r-', label='C-cell')
plt.plot(t, sol[:, 3], 'm-', label='M.A.-cell')
plt.plot(t, sol[:, 4], 'b--', label='mRNA')
plt.legend(loc='best')
plt.xlabel('t(s)')
plt.ylabel('C(mM)')
plt.grid()
plt.show()


# In[21]:


plt.plot(t, sol[:, 5], 'g--', label='LigM')
plt.plot(t, sol[:, 6], 'r--', label='AroY')
plt.plot(t, sol[:, 7], 'm--', label='CatA')
plt.legend(loc='best')
plt.xlabel('t(s)')
plt.ylabel('C(mM)')
plt.grid()
plt.show()


# In[22]:


print("Final concentrations of:")
print("ligM subunit: "+str(round(sol[400003, 5],5))+" mM")
print("ligM enzyme: "+str(round((sol[400003,5])/2,5))+" mM")
print("aroY subunit: "+str(round(sol[400003, 6],5))+" mM")
print("aroY enzyme: "+str(round((sol[400003,6])/6,5))+" mM")
print("catA subunit: "+str(round(sol[400003, 7],5))+" mM")
print("catA enzyme: "+str(round((sol[400003,7])/2,5))+" mM")


# In[23]:


print("Final mass of:")
print("ligM: "+str(round((sol[400003,5])*55325/(1000),1))+" fg")
print("aroY: "+str(round((sol[400003,6])*54012/(1000),1))+" fg")
print("catA: "+str(round((sol[400003,7])*33434/(1000),1))+" fg")
print()
print("Common E. Coli dry weight: 284 fg")
print("Common total protein mass on E. Coli: 156 fg")
print("Most expressed protein mass on common E. Coli: ~7.8 fg (MetE)")


# In[ ]:




