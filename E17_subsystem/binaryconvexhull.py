# In[1]
# 导入模块
#import matplotlib.pyplot as plt
import matplotlib  # 解决负数显示问题
from ase.phasediagram import PhaseDiagram  # 用于绘制相图
from getConvexHull import read_convexhull  # 用于读取凸壳数据
# In[2]
matplotlib.rcParams['axes.unicode_minus']=False # 解决负数显示问题
# In[3]
filepath_0 = 'E17_subsystem/input_files/extended_convex_hull_0GPa_prepared' # 0GPa凸壳数据文件路径

subsystem_LiN = [0,2] # LiN子系统
subsystem_LiC = [0,1] # LiC子系统
subsystem_CN = [1,2] # CN子系统

datastart = 7 # 数据起始行
filetype = 'extended_convex_hull' # 文件类型

elements = ['Li','C','N'] # 元素列表

atomenergy_0 = [-1.904144347, -9.225171725, -8.320979825] # 原子能量
atomenergy_50 = [1.710545228125, -7.40000050875, -4.917002715]
atomenergy_100 = [3.909955566, -5.843035403, -3.063771986]

fitreq = 0 # fitness要求

ref_E = [
        ('Li', 0),
        #('C', 0),
        ('N', 0)
         ]
# In[4]
# 读取凸壳数据
convexhull = read_convexhull(filepath=filepath_0,
                                 filetype=filetype,
                                 elements=elements,
                                 datastart=datastart,
                                 atomenergy=atomenergy_0,
                                 subsystem=subsystem_LiN,
                                 fitreq=fitreq
                                 )
# In[5]
# 加载数据
convexhull.loadfile()
convexhull.selectByFitness() # 选择fitness满足要求的数据
# In[6]
convexhull.selectSubsystem() # 选择子系统
convexhull.getToten() # 获取总能
convexhull.getHf() # 获取形成能
refs = convexhull.getRefs() # 用于绘制相图
#print(refs)

ref_LiN_mp_0=[('Li1N3', -0.146 * 4),
              ('Li3N1', -0.399 * 4)]
    
refs_uspex = refs
refs_uspex.extend(ref_E)
refs_uspex.extend(ref_LiN_mp_0)
# In[7]
pd = PhaseDiagram(refs_uspex)
pd.plot(show=True,dims=2)

# ref_LiN_article_0 = [('Li8N8', -5.709), ('Li2N2', -1.4515), ('Li8N8', -5.292),
#                      ('Li3N1', -1.6214), ('Li2N4', -1.498), ('Li2N6', -1.02)]
# ref_LiN_article_50 = [('Li2N2', -6.190), ('Li8N8', -26.7942), ('Li12N8', -34.62690302), ('Li6N4', -17.2129),
#                       ('Li3N1', -2.2998), ('Li26N2', -14.6908), ('Li2N4', -8.0564), ('Li2N6', -6.9579), ('Li4N20', -20.4620)]
# ref_LiN_article_100 = [('Li8N8', -32.3483), ('Li12N8', -42.7455), ('Li5N1', -9.337723944), ('Li3N1', -0.7045),
#                        ('Li26N2', -18.1087), ('Li2N4', -9.4224), ('Li2N6', -5.9945), ('Li4N20', -20.3890)]
# ref_LiC_article_0 = [('Li4C4', 0.0054), ('Li12C6', -0.0170),
#                      ('Li8C6', -0.5732), ('Li1C6', -0.0885), ('Li2C24', -0.2381)]
# ref_LiC_article_50 = [('Li12C6', -9.396), ('Li8C6', -6.124), ('Li16C16', -15.186), ('Li12C12', -11.3251), ('Li8C12', -9.2844), ('Li16C8', -15.0334),
#                       ('Li6C8', -6.8365), ('Li8C2', -4.7572), ('Li18C3', -7.4144), ('Li24C9', -19.5666), ('Li32C4', -10.1814), ('Li2C12', 0.52)]
# ref_LiC_article_100 = [('Li12C6', -13.5896), ('Li8C6', -6.4144), ('Li16C16', -20.16), ('Li12C12', -14.5725), ('Li8C12', -12.2047),
#                        ('Li16C8', -20.256), ('Li6C8', -8.805), ('Li8C2', -7.5857), ('Li18C3', -11.1642), ('Li24C9', -29.962), ('Li32C4', -15.464), ('Li2C12', 2.035)]
# ref_CN_article_0 = []
# ref_CN_article_50 = [('C4N4', -2.243), ('C4N4', -2.715),
#                      ('C8N16', -7.195), ('C12N16', -11.11), ('C6N8', -5.129)]
# ref_CN_article_100 = [('C4N4', -2.668), ('C4N4', -2.717),
#                       ('C8N16', -10.26), ('C12N16', -12.51), ('C6N8', -6.725)]
# ref_referance = ref_LiN_article_0
# ref_referance.extend(ref_E)
# pd = PhaseDiagram(ref_referance)
# pd.plot(show=True, dims=2)
# refs.extend(ref_referance)
# pd = PhaseDiagram(refs)
# pd.plot(show=True, dims=2)
