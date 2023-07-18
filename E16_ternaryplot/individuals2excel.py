# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 12:55:07 2023

@author: 27682
"""

# In[1]
from getConvexHull import read_convexhull
import pandas as pd

# In[2]
# basic input parameters
filepath = r'C:/Users/27682/PycharmProjects/uspex_tools/E16_ternaryplot/input/Individuals_prepared' # using absolute path here
datastart = 3
filetype = 'individuals'
atomenergy = [1.710545228125, -7.40000050875, -4.917002715]  # LiCN 50GPa
elements = ['Li', 'C', 'N']

# In[3]
# load file
convexhull = read_convexhull(filepath=filepath,
                             filetype=filetype,
                             elements=elements,
                             atomenergy=atomenergy,
                             datastart=datastart)
convexhull.loadfile()

# In[4]
#calculate energy
convexhull.getToten()
convexhull.getHf()

# In[5]
# remove duplicate
convexhull.remove_dup()

# In[6]
# convert your data to excel

# composition normalize
comp_normed_vector = []
for comp in convexhull.compositions_vector:
    comp_normed_vector.append([x / sum(comp) for x in comp])
# to separate compositions into different columns
x = []
y = []
z = []
for comp in comp_normed_vector:
    x.append(comp[1])
    y.append(comp[2])
    z.append(comp[0])
df = pd.DataFrame({'comp':convexhull.compositions,'C':x,'N':y,'Li':z,"hf":convexhull.Hf})
df.to_excel(r'C:/Users/27682/PycharmProjects/uspex_tools/E16_ternaryplot/output/Individuals_50.xlsx')