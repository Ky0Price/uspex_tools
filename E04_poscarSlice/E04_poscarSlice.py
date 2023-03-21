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
# 注意：
# 1. 此处的文件路径是相对于本脚本的路径，而不是相对于python的路径
# 2. 此处用的是经过一定预处理后的extended_convex_hull文件,可以看见文件名为extended_convex_hull_0GPa_prepared
#   该文件是在extended_convex_hull（即extended_convex_hull_0GPa_original）文件的基础上，
#   将一些numpy无法识别的字符替换为了空格，以便于读取
filetype = 'extended_convex_hull'
elements = ['Li','C','N']
pospath = 'input_files/extended_convex_hull_POSCARS_0GPa'
# 存放总POSCARS集合文件的路径
pos_savepath = 'poscars'
# 存放单个POSCAR文件的目标路径
potpath = 'input_files'
# 存放各种POTCAR文件的路径
fitreq = 0.0

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
# 开始切割POSCAR
convexhull.poscar_slice()
# 根据fitreq，筛选目标结构
convexhull.poscar_selectByFitness()
# 生成POSCAR文件，以文件夹形式存放在pos_savepath中
convexhull.poscar_buildfile()

# In[6]
# 放入对应的POTCAR文件
convexhull.getPOTCAR()

# In[7]
# 生成INCAR文件
# 未来也许会添加这个功能，但目前感觉就bash脚本一行代码的事情，没必要写进来
