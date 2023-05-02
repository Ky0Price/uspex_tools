'''
#!~/bin/python
#-*- coding: UTF-8 -*-
'''
def FormationEnCalc(idlist,system,moleEn,atomEn,op,pattern):
    # idlist对应rC中getID()返回的结构ID列表,没发现有什么用，但是删掉会出问题，就给他留着吧
    # moleEn指分子总焓，对应的是rC中getEnthalpies()返回的总焓列表
    # system指体系的组成向量，对应rC中getTotalEnt()返回的组成向量列表
    # atomEn指的是组成分子前每个原子的基态能量，需要自己手动写个字典传入，比如{Li:-2,C:-9}
    # 形成焓计算公式deltaE(mAnB) = E(mAnB)-m*E(A)-n*E(B)
    FElist=[] #存放最终的ID-生成焓的对应字典的列表
    if pattern == 'scan':
        for id,comp,en in zip(idlist,system,moleEn):
            # 计算体系中原子种类数和各原子种类对应的原子数目
            FEdict = {}
            FormationEn = en
            complist = comp.split(' ') # 按元素切割一下
            #print(complist)
            totalAtomNumber = 0
            composition = ''  # 用于生成方便输入ASE的组成格式
            for ele in complist[:-1]:
                element = ''  # 记录元素种类
                number = ''  # 记录原子个数
                #print(ele)
                for i in ele:
                    if i.isalpha():
                        element += i
                    else:
                        number += i
                #print(element)
                #print(number)
                totalAtomNumber += eval(number)
                composition += element+number
                FormationEn -= atomEn[element]*eval(number) # 减去对应原子的基态能量
            #print(composition)
            print(totalAtomNumber)
            if op==1:
                FEdict['composition']=composition # ASE要求输入的是分子总生成焓，如果想要得到分子中每原子的平均生成焓，请把这一项改为“FEdict[composition]=FormationEn/totalAtomNumber”
                FEdict['FE']=FormationEn
                FEdict['id']=id
                #print(FEdict)
            else:
                FEdict['composition'] = composition
                FEdict['FE'] = FormationEn/totalAtomNumber
                FEdict['id'] = id
            FElist.append(FEdict)
    else:
        for id,comp,en in zip(idlist,system,moleEn):
            # 计算体系中原子种类数和各原子种类对应的原子数目
            FEdict = {}
            FormationEn = en
            complist = comp.split(' ') # 按元素切割一下
            #print(complist)
            totalAtomNumber = 0
            composition = ''  # 用于生成方便输入ASE的组成格式
            for ele in complist:
                element = ''  # 记录元素种类
                number = ''  # 记录原子个数
                #print(ele)
                for i in ele:
                    if i.isalpha():
                        element += i
                    else:
                        number += i
                #print(element)
                #print(number)
                totalAtomNumber += eval(number)
                composition += element+number
                FormationEn -= atomEn[element]*eval(number) # 减去对应原子的基态能量
            #print(composition)
            if op==1:
                FEdict['composition']=composition # ASE要求输入的是分子总生成焓，如果想要得到分子中每原子的平均生成焓，请把这一项改为“FEdict[composition]=FormationEn/totalAtomNumber”
                FEdict['FE']=FormationEn
                FEdict['id']=id
                #print(FEdict)
            else:
                FEdict['composition'] = composition
                FEdict['FE'] = FormationEn/totalAtomNumber
                FEdict['id'] = id
            FElist.append(FEdict)
    return FElist

def FE2refs(FElist):
    refs=[]
    for entry in FElist:
        ref=(entry['composition'],entry['FE'])
        refs.append(ref)
    return refs
    

def getSpaceGroup(path,symtole=1e-1):
    # path指存放POSCAR的文件夹目录的列表
    # 本脚本暂时无法处理子文件夹
    # symtole指找对称群时的tolerance
    import os

    import ase.io.vasp
    import ase.spacegroup

    #遍历文件夹中的文件
    sglist = []
    for i in path:
        for root,dirs,files in os.walk(i):
        #获取ID
        #print(files)
        #记录走入的文件夹层数
            for file in files:
                #print(file)
                if file[0] == 'E':
                    id=''
                    for i in file:
                        if i == '_':
                            break
                        id += i
                elif file.split('_')[1][0] == 'm':
                    id=''
                    for i in file.split('_'):
                        if i[0] == 'm':
                            id += i
                else:
                    id = file.split('_')[0]
                file = root+'/'+file #生成完整路径
                #print(file)
                atom = ase.io.vasp.read_vasp(file)
                #print(atom)
                sg = ase.spacegroup.get_spacegroup(atom,symprec=symtole)
                #print(sg)
                sgd=sg.todict()
                #print(sgd)
                sgd['id']=id
                sglist.append(sgd)
    #print(len(sglist))
    #print(sglist)
    return sglist


    
def datawirting_exe(datapath,comp_FE_id_list=[],pressure=0,spacegrouplist=[],Enthalpylist=[]):
    # 整理数据用
    # datapath写入你希望保存文件的位置
    # idlist传入结构的id，如EA123或者mp-1234或123456都行
    # compAndFEdict传入最后的ref（refs，refs），他包含了结构和对应的分子形成焓
    # spacegrouplist指结构的空间群
    # pressure默认0Gpa
    # enthalpylist传入分子的总焓
    import csv
    with open(datapath,'w',encoding='utf-8',newline='') as f:
        writer = csv.writer(f)
        # 标题行
        head = ['id',
                'compsition',
                'spacegroup',
                'pressure(GPa)',
                'enthalpy    (ev/atom)',
                'formation enthalpy    (ev/atom)']
        writer.writerow(head)
        for comp_FE_id,enthalpy in zip(comp_FE_id_list,Enthalpylist):
            comp = comp_FE_id['composition']
            FE = comp_FE_id['FE']
            id = comp_FE_id['id']
            #print(id)
            pre = pressure
            sg = 0
            for sgdict in spacegrouplist:
                if sgdict['id'] == id:
                    #print(sgdict['number'])
                    sg += sgdict['number']
            #print(sg)
            row = [id]+[comp]+[sg]+[pre]+[enthalpy]+[FE]
            writer.writerow(row)


          
import pandas as pd
from ase.phasediagram import PhaseDiagram

data_file='./data.csv'
df = pd.read_csv(data_file,names=['comp','energy'])
#print(df['comp'])
idlist=[]
system=[]
for i in df['comp']:
    idlist.append(i.split('_')[0])
    system.append(' '.join(i.split('_')[1:]))

moleEn=df['energy']
#print(moleEn)

Li = -2.342349913

C = -10.09022995

N = -9.288053935

atomEn = {'Li':Li,'C':C,'N':N}
#print(idlist)
#print(system)
FElist_S = FormationEnCalc(idlist=idlist,moleEn=moleEn,system=system,atomEn=atomEn,op=2,pattern='scan')
refs=FE2refs(FElist=FElist_S)
ref_E=[('Li',0),('C',0),('N',0)]
refs.extend(ref_E)
pd = PhaseDiagram(refs)
pd.plot(show=True,dims=2)


idlist_M=['mp-1001079','mp-30057','mp-1190169','mp-1190940','mp-9610','mp-1247387','mp-2251','mp-1021323','mp-2659']
idlist_M=['mp-1001079','mp-30057','mp-1190169','mp-1190940','mp-9610','mp-1247387','mp-2251','mp-1021323','mp-2659']
idlist_Z = ['mp-1018134','mp-569304','mp-25']
mpsys=['Li4 C8 N8','Li4 C4 N4','Li2 C12 N8','Li4 C8 N12','Li4 C2 N4','Li12 C4 N8','Li3 N1','Li1 C12','Li2 N6']
#print(mpsystem) #输出没问题
mpmoleEn=[-162.63248909,-91.15394848,-194.98214572,-205.23673597,-74.77445857,-155.69452599,-18.27389958,-123.78533053,-62.09848495]
mpmoleEn_a=[-162.63248909/20,-91.15394848/12,-194.98214572/22,-205.23673597/24,-74.77445857/10,-155.69452599/24,-18.27389958/4,-123.78533053/13,-62.09848495/8]
FElist_M=FormationEnCalc(idlist=idlist_M,moleEn=mpmoleEn,system=mpsys,atomEn=atomEn,op=2,pattern='mp')
refMP=FE2refs(FElist_M)
#print(refMP)
#refMP.extend(ref_E)
refs.extend(refMP)
#pd = PhaseDiagram(refs)
#pd.plot(show=True,dims=2)

idlist=idlist+idlist_Z+idlist_M
FElist_Z=[{'composition':'Li','FE':0,'id':'mp-1018134'},
          {'composition':'C','FE':0,'id':'mp-569304'},
          {'composition':'N','FE':0,'id':'mp-25'}]
comp_FE_id_list = FElist_S+FElist_Z+FElist_M
Enthalpylist=[]
for i in moleEn:
    Enthalpylist.append(i)
for i in atomEn.values():
    Enthalpylist.append(i)
Enthalpylist_M = mpmoleEn_a
Enthalpylist += Enthalpylist_M
uspexpath='./uspex_pos'
mppath='./mp_pos'
elepath='./ele_pos'
dirlist=[uspexpath,elepath,mppath]
pressure=0
datapath=['./scandata_0GPa_sym1e-1.csv']
symtole_list=[1e-1]
for i,j in zip(datapath,symtole_list):
    spacegruoplist = getSpaceGroup(path=dirlist, symtole=j)
    datawirting_exe(datapath=i,comp_FE_id_list=comp_FE_id_list,pressure=pressure,spacegrouplist=spacegruoplist,Enthalpylist=Enthalpylist)


    
    

