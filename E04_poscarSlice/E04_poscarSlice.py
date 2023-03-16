# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:34:24 2023

@author: 27682
"""

# In[1]
#导入依赖库
from getConvexHull import read_convexhull

# In[2]
#主要输入参数准备,以下路径均为相对于本脚本的路径
filepath = 'input_files/extended_convex_hull_0GPa_prepared'
filetype = 'extended_convex_hull'
elements = ['Li','C','N']
pospath = 'input_files/extended_convex_hull_POSCARS_0GPa'
# 存放总POSCARS集合文件的路径
pos_savepath = 'poscars'
# 存放单个POSCAR文件的目标路径
potpath = 'input_files'
# 存放各种POTCAR文件的路径
fitreq = 0.03

# In[3]
#创建对象
convexhull = read_convexhull(
                             filepath = filepath, 
                             filetype = filetype, 
                             elements = elements,
                             pospath = pospath,
                             pos_savepath = pos_savepath,
                             potpath=potpath,
                             fitreq = fitreq
                             )

# In[4]
# 根据提供的fitreq，筛选目标结构
convexhull.loadfile()
convexhull.selectByFitness()

# In[5]
# 开始切割POSCAR，筛选，最后创建单独文件夹
convexhull.poscar_slice()
convexhull.poscar_selectByFitness()
convexhull.poscar_buildfile()

# In[6]
# 放入对应的POTCAR文件
convexhull.getPOTCAR()

