#!/usr/bin/env python
# coding: utf-8

#Metabolic pathwey model with modified BMC
# In[1]:


def pend(z, t, K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, j_m_V, j_m_P, j_m_C, j_m_M, j_b_Vo, j_b_Po, j_b_Co, j_b_Mo, j_b_Vi, j_b_Pi, j_b_Ci, j_b_Mi, Vo, Po, Co, Mo):
    Vc, Pc, Cc, Mc, Vb, Pb, Cb, Mb = z #Variables to integrate
    dVcdt = (j_b_Vi*(Vb-Vc)+j_b_Vo*(Vb-Vc)-j_m_V*(Vc-Vo))/volc-(K1cat*Ec*Vc)/(Vc + K1m) #Vanilic acid
    dPcdt = (j_b_Pi*(Pb-Pc)+j_b_Po*(Pb-Pc)-j_m_P*(Pc-Po))/volc+(K1cat*Ec*Vc)/(Vc + K1m) - (K2cat*EEc*Pc)/(Pc + K2m) #PCA
    dCcdt = (j_b_Ci*(Cb-Cc)+j_b_Co*(Cb-Cc)-j_m_C*(Cc-Co))/volc+(K2cat*EEc*Pc)/(Pc + K2m) - (K1cat*EEEc*Cc)/(Cc + K3m) #Catechol
    dMcdt = (j_b_Mi*(Mb-Mc)+j_b_Mo*(Mb-Mc)-j_m_M*(Mc-Mo))/volc+(K1cat*EEEc*Cc)/(Cc + K3m) #Cis,Cis-Muconic Acid
    dVbdt = (-j_b_Vi*(Vb-Vc)-j_b_Vo*(Vb-Vc))/volb #Vanilic acid in BMC
    dPbdt = (-j_b_Pi*(Pb-Pc)-j_b_Po*(Pb-Pc))/volb - (K2cat*EEb*Pb)/(Pb + K2m) #PCA in BMC
    dCbdt = (-j_b_Ci*(Cb-Cc)-j_b_Co*(Cb-Cc))/volb+(K2cat*EEb*Pb)/(Pb + K2m) - (K3cat*EEEb*Cb)/(Cb + K3m) #Catechol in BMC
    dMbdt = (-j_b_Mi*(Mb-Mc)-j_b_Mo*(Mb-Mc))/volb+(K3cat*EEEb*Cb)/(Cb + K3m) #Muconic acid in BMC
    dCdt = [dVcdt,dPcdt,dCcdt,dMcdt,dVbdt,dPbdt,dCbdt,dMbdt] #Functions to integrate
    return dCdt


# In[2]:


K1cat = (0.078*9.14) #Kcat LigM
K1m = 0.078 #Km LigM
K2cat = 0.00009*1 #Kcat AroY
K2m = 0.59 #Km AroY
K3cat = (0.00185*10485) #Kcat CatA
K3m = 0.00185 #Kmm CatA
volo = 1000000 #Medium volume (outside)
volc = 1e-15 #Citoplasm volume
volb = 0.1e-15 #BMC volume
Vo = 10e-3 #vanilic acid outside
Po = 0 #protocatechuic acid outside
Co = 0 #catechol outside
Mo = 0 #cis,cis-muconic acid outside
Ec = 0.1479
EEc = 0.04329
EEEc = 0.2039
Eb = 1.2*Ec
EEb = 1.2*EEc
EEEb = 1.2*EEEc
#flux constants for species trough membrane
j_m_V = 1e-13 # VA trough Membrane
j_m_P = 1e-13 # PCA
j_m_C = 1e-13 # Catechol
j_m_M = 1e-13 # MA
j_b_Vo = 0.1e-13 #through pduU
j_b_Po = 0.1e-13 
j_b_Co = 0.1e-13
j_b_Mo = 0.1e-13
j_b_Vi = 1e-15 #pduA
j_b_Pi = 1e-15
j_b_Ci = 1e-15
j_b_Mi = 1e-15


# In[3]:


z0 = [0,0,0,0,0,0,0,0] #Initial values


# In[4]:


import numpy as np
t = np.linspace(0, 0.1, 100001) #Create 100001 values between 0-0.1 in the X axis


# In[5]:


from scipy.integrate import odeint #Integrate all functions for the X axis given values
sol1 = odeint(pend, z0, t, args=(K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, j_m_V, j_m_P, j_m_C, j_m_M, j_b_Vo, j_b_Po, j_b_Co, j_b_Mo, j_b_Vi, j_b_Pi, j_b_Ci, j_b_Mi, Vo, Po, Co, Mo))


# In[6]:


z0 = list(sol1[100000,:]) #New initial values using final values for past function


# In[7]:


t = np.linspace(0.1, 10, 100001) #Give new values to the X axis for the next step


# In[8]:


from scipy.integrate import odeint
sol2 = odeint(pend, z0, t, args=(K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, j_m_V, j_m_P, j_m_C, j_m_M, j_b_Vo, j_b_Po, j_b_Co, j_b_Mo, j_b_Vi, j_b_Pi, j_b_Ci, j_b_Mi, Vo, Po, Co, Mo))


# In[9]:


z0 = list(sol2[100000,:])


# In[10]:


t = np.linspace(10, 1000, 100001)


# In[11]:


from scipy.integrate import odeint
sol3 = odeint(pend, z0, t, args=(K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, j_m_V, j_m_P, j_m_C, j_m_M, j_b_Vo, j_b_Po, j_b_Co, j_b_Mo, j_b_Vi, j_b_Pi, j_b_Ci, j_b_Mi, Vo, Po, Co, Mo))


# In[12]:


z0 = list(sol3[100000,:])


# In[13]:


t = np.linspace(1000, 1800, 100001)


# In[14]:


from scipy.integrate import odeint
sol4 = odeint(pend, z0, t, args=(K1cat, K1m, K2cat, K2m, K3cat, K3m, volo, volc, volb, j_m_V, j_m_P, j_m_C, j_m_M, j_b_Vo, j_b_Po, j_b_Co, j_b_Mo, j_b_Vi, j_b_Pi, j_b_Ci, j_b_Mi, Vo, Po, Co, Mo))


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
d = np.linspace(1000, 1800, 100001)
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
plt.legend(loc='best')
plt.xlabel('t(s)')
plt.ylabel('C(mM)')
plt.grid()
plt.show()


# In[20]:


plt.plot(t, sol[:, 1], 'g-', label='PCA-cell')
plt.plot(t, sol[:, 2], 'r-', label='C-cell')
plt.plot(t, sol[:, 3], 'm-', label='M.A.-cell')
plt.legend(loc='best')
plt.xlabel('t(s)')
plt.ylabel('C(mM)')
plt.grid()
plt.show()


# In[21]:


plt.plot(t, sol[:, 4], 'b--', label='V.A.-BMC')
plt.plot(t, sol[:, 5], 'g--', label='PCA-BMC')
plt.plot(t, sol[:, 6], 'r--', label='C-BMC')
plt.plot(t, sol[:, 7], 'm--', label='M.A.-BMC')
plt.legend(loc='best')
plt.xlabel('t(s)')
plt.ylabel('C(mM)')
plt.grid()
plt.show()


# In[22]:


print("Final concentrations of:")
print("Vanilic acid:")
print("Cell- "+str(round(sol[400003, 0],5))+"mM")
print("BMC- "+str(round(sol[400003, 4],5))+"mM")
print()
print("PCA:")
print("Cell- "+str(round(sol[400003, 1],6))+"mM")
print("BMC- "+str(round(sol[400003, 5],6))+"mM")
print()
print("Catechol:")
print("Cell- "+str(round(sol[400003, 2],12))+"mM")
print("BMC- "+str(round(sol[400003, 6],12))+"mM")
print()
print("Muconic acid:")
print("Cell- "+str(round(sol[400003, 3],14))+"mM")
print("BMC- "+str(round(sol[400003, 7],13))+"mM")
print()


# In[ ]:




