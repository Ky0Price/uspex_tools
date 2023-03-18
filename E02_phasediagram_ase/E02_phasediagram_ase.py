# 处理USPEX输出文件，并利用ase库绘制相图
# 作者：zzn
# 示例脚本，分cell代码块


from os import mkdir, path  # 用于创建文件夹

import matplotlib  # 绘图库
from ase import phasediagram  # 用于绘制相图
from FormationEnCalc import FE2refs, FormationEnCalc  # 计算外部数据的形成焓
# In[1]: 导入依赖库
from getConvexHull import read_convexhull  # 读取USPEX输出文件

# In[2]: 解决乱码问题
matplotlib.rcParams['font.sans-serif'] = ['SimHei'] # 指定默认字体
matplotlib.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号'-'显示为方块的问题

# In[3]: 读取USPEX输出文件
# 读取USPEX输出文件，得到相图数据，这里以Li-C-N体系为例
filepath = "E02_phasediagram_ase\input_files\extended_convex_hull_0GPa_prepared"
# 注意：
# 1. 此处的文件路径是相对于本脚本的路径，而不是相对于python的路径
# 2. 此处用的是经过一定预处理后的extended_convex_hull文件,可以看见文件名为extended_convex_hull_0GPa_prepared
#   该文件是在extended_convex_hull（即extended_convex_hull_0GPa_original）文件的基础上，
#   将一些numpy无法识别的字符替换为了空格，以便于读取
filetype = 'extended_convex_hull'
elements = ['Li', 'C', 'N']
atom_energy = [-1.904144347, # Li单质在0GPa下的能量，eV/atom
                -9.225171725, # C单质在0GPa下的能量，eV/atom
                -8.320979825] # N单质在0GPa下的能量，eV/atom
                # 各元素单质的能量/焓值，
                # 用于计算形成焓（详情参考）E00的第二小节，
                # 此处能量值的顺序要与elements中元素的顺序一一对应；
                # 由于本次示例中extended_convex_hull来自0GPa下的结果，
                # 因此此处要使用0GPa下各元素单质的能量
datastart = 7
fitreq = 0
ref_E = [('Li',0),('C',0),('N',0)]

# In[4]: 读取数据
# 创建对象
convexhull = read_convexhull(
                            filepath=filepath, 
                            filetype=filetype, 
                            elements=elements, 
                            datastart=datastart, 
                            fitreq=fitreq)
# 加载数据
convexhull.loadfile()

# In[5]: 绘制相图
# 根据fitreq参数，筛选数据
convexhull.selectByFitness()
# 获得能量数据
convexhull.getToten()
convexhull.getHf()
# 获得用于ase绘图的数据
refs_upsex = convexhull.getRefs()
# 加入单质
refs_upsex.extend(ref_E)

# In[6]: 绘制相图
# 创建绘图对象
pd = phasediagram.PhaseDiagram(refs_upsex)
# 指定相图维度
dims = 2 # 2D相图
# dims = 3 # 3D相图
pd.plot(dims = dims)
# 保存相图
while not path.exists('./output_files'):
    try:
        mkdir('./output_files')
    except:
        pass
pd.savefig('./output_files/phase_diagram.png') 

# In[7]: 放入参考数据
# 放入从数据库或者参考文献中获得的参考数据，用于验证uspex计算结果的准确性和完备性
# 这里以materials project数据库中的数据为例
# 创建id列表，用于辨识不同的化合物
mp_idlist = [
            'mp-1001079', 'mp-30057', 'mp-1190169', 'mp-1190940',
            'mp-9610', 'mp-1247387', 'mp-2251', 'mp-1021323', 'mp-2659'
            ]
# 创建化合物列表，用于存放化合物对象
mp_compounds = [
                'Li4_C8_N8', 'Li4_C4_N4', 'Li2_C12_N8', 'Li4_C8_N12',
                'Li2_C1_N2', 'Li12_C4_N8', 'Li3_N1', 'Li1_C12', 'Li2_N6'
                ]
# 创建能量列表，用于存放化合物的能量
mp_toten = [
            -145.36960, -81.086877, -175.80325, -183.56898, -6.638*5, 
            -136.75012, -15.627720, -112.7152, -54.901334
            ]
# 各个元素单质的能量，以字典形式存储（代码的年代原因，数值与上面atom_energy一致，但是格式不一样，日后会改进）
mp_atom_energy = {
                'Li': -1.904144347,
                'C': -9.225171725,
                'N': -8.320979825
}

# In[8]: 计算参考数据的形成焓
ref_mp = FE2refs(
    FormationEnCalc(
        idlist=mp_idlist,
        system=mp_compounds,
        moleEn=mp_toten,
        atomEn=mp_atom_energy
    )
)

# In[9]: 与uspex数据融合
# 将uspex数据与参考数据融合
refs = refs_upsex + ref_mp

# In[10]: 绘制相图
# 创建绘图对象
pd = phasediagram.PhaseDiagram(refs)
# 指定相图维度
dims = 2 # 2D相图
# dims = 3 # 3D相图
pd.plot(dims = dims)
# 保存相图
pd.savefig('./output_files/phase_diagram_with_refs.png')
