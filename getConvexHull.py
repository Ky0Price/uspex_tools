# -*- coding: utf-8 -*-
"""
Created on Fri Oct 1 10:05:20 2021

@author: 张圳南@湖南大学
"""
class read_convexhull:
    # 该类主要用于USPEX输出文件的后处理，包括extended_convex_hull和individuals文件
    # 目前能进行的功能有：
    # 1. 读取extended_convex_hull文件
    # 2. 读取individuals文件
    # 3. 读取poscar总文件并切片输出
    # 4. 计算总能量，形成焓，Gibbs能
    # 5. 生成ase画相图所需要的数据格式refs
    # 6. 去掉fitness大于某个值的结构
    # 7. 选取子系统的结构
    # 8. 去掉化学式相同但是能量更高的结构

    
    def __init__(self,filepath,filetype,elements,datastart=7,pospath=None,pos_savepath=None,potpath=None,atomenergy=[],subsystem=None,fitreq=0.0):
        # 默认extended_convex_hull第7行开始正式有数据
        # 行数请根据实际文件进行设置
        self.filepath = filepath  # 数据文件的目录
        self.pospath = pospath
        self.filetype = filetype
        self.fitreq = fitreq  # 获取的fitness的最大值
        self.datastart = datastart  # 数据开始的行数
        self.elements = elements
        self.elenum = len(elements)
        self.subsystem = subsystem
        self.atomenergy = atomenergy
        self.pos_savepath = pos_savepath
        self.potpath = potpath
        #self.posreq = posreq
        self.data = []
        self.id = []
        self.id_sub = []
        self.compositions = []
        self.compositions_vector = []
        self.compositions_sub = []
        self.compositions_vector_sub = []
        self.subsystem_loc = []
        self.enthalpy = []
        self.enthalpy_sub = []
        self.fitness = []
        self.fitness_sub = []
        self.volume = []
        self.selectlen = 0
        self.sym = []
        self.totalEnergy = []
        self.totalEnergy_sub = []
        #self.totalEnergy_per_atom = []
        self.Hf = []
        #self.Hf_per_atom = []
        self.refs = []
        self.poscar_silced = []
        self.poscar_silced_selected = []
        self.poscarpaths = []
        self.gibbs = []
        
        



    
    def loadfile(self):
        # 读取文件，将数据存入对象的基本属性
        # 对象的基本属性
        # ：id、组成、总焓、体积、fitness
        import numpy as np
        np.set_printoptions(precision=None,suppress=True)
        self.data=np.loadtxt(self.filepath,skiprows=self.datastart-1)
        if self.filetype=='extended_convex_hull':
            for i in self.data:
                self.id.append(int(i[0]))
            entry_num=0
            for entry in self.data:
                entry_num+=1
                composition=''
                composition_vector=[]
                for i in range(self.elenum):
                    composition+=self.elements[i]+str(int(entry[i+1]))
                    composition_vector.append(int(entry[i+1]))
                self.compositions.append(composition)
                self.compositions_vector.append(composition_vector)
            self.enthalpy=self.data[:,self.elenum+1]
            self.volume=self.data[:,self.elenum+2]
            self.fitness=self.data[:,self.elenum+3]
        elif self.filetype=='individuals':
            self.id=self.data[:,1]
            for entry in self.data:
                composition=''
                composition_vector = []
                for i in range(self.elenum):
                    composition+=self.elements[i]+str(int(entry[i+2]))
                    composition_vector.append(int(entry[i+2]))
                self.compositions.append(composition)
                self.compositions_vector.append(composition_vector)
            self.enthalpy=self.data[:,self.elenum+2]
            self.volume=self.data[:,self.elenum+3]
            self.fitness=self.data[:,self.elenum+4]
        return  self

    def selectByFitness(self):
        # 依据fitness的值选取结构，上限为用户输入的fitreq
        # 感觉这里最好是新建个对象，而不是直接修改原来的对象，也就是弄个子类
        # 但是我懒得弄了，就这样吧，反正能用
        if self.filetype!='extended_convex_hull':
            print('ERR! To use this function, filetype should be extended_convex_hull')
            return
        else:
            for i in self.fitness:
                while i <= self.fitreq:
                    self.selectlen+=1
                    break
        self.id=self.id[:self.selectlen]
        self.compositions=self.compositions[:self.selectlen]
        self.compositions_vector=self.compositions_vector[:self.selectlen]
        self.enthalpy=self.enthalpy[:self.selectlen]
        return self

    def selectSubsystem(self):
    # 获取子系统的结构
    # 什么是子系统？
    # 子系统是指用户自己定义的一部分原子，
    # 比如用户想要从一个C-H-O-N的系统中获取含有C、H、O的结构，那么这个子系统就是C、H、O
    # 在这里，我们利用subsystem这个参数来定义子系统
    # subsystem是一个列表，列表中的元素是element中对应元素的序号
        compnum=-1
        if self.filetype=='extended_convex_hull':
            for comp in self.compositions_vector:
                compnum+=1
                composition_sub=''
                composition_vector_sub=[]
                jud = 1 - 1 * len(self.elements)
                for i in range(len(self.elements)):
                    if comp[i] != 0 and i in self.subsystem:
                        jud += 1
                    elif comp[i] == 0 and i not in self.subsystem:
                        jud += 1
                    else:
                        jud += 0
                if jud == 1:
                    for i in self.subsystem:
                        composition_sub += self.elements[i] + str(int(comp[i]))
                        composition_vector_sub.append(int(comp[i]))
                    self.compositions_sub.append(composition_sub)
                    self.compositions_vector_sub.append(composition_vector_sub)
                    self.subsystem_loc.append(compnum)
                    self.enthalpy_sub.append(self.enthalpy[compnum])
                    self.fitness_sub.append(self.fitness[compnum])
                    self.id_sub.append(self.id[compnum])
        elif self.filetype=='individuals':
            for comp in self.compositions_vector:
                compnum += 1
                composition_sub = ''
                composition_vector_sub = []
                jud = 1 - 1 * len(self.elements)
                for i in range(len(self.elements)):
                    if comp[i] != 0 and i in self.subsystem:
                        jud += 1
                    elif comp[i] == 0 and i not in self.subsystem:
                        jud += 1
                    else:
                        jud += 0
                if jud == 1:
                    for i in self.subsystem:
                        composition_sub += self.elements[i] + str(int(comp[i ]))
                        composition_vector_sub.append(int(comp[i]))
                    self.compositions_sub.append(composition_sub)
                    self.compositions_vector_sub.append(composition_vector_sub)
                    self.subsystem_loc.append(compnum)
                    self.enthalpy_sub.append(self.enthalpy[compnum])
                    self.fitness_sub.append(self.fitness[compnum])
                    self.id_sub.append(self.id[compnum])
        return self


    def getToten(self):  
        # 如果文件形式是extended_convex_hull，
        # 输出以化学式为单位的整体能量；
        # 如果是individuals，则输出单个原子的平均能量
        if self.subsystem==None:
            for en,comp in zip(self.enthalpy,self.compositions_vector):
                totatom=sum(comp)
                if self.filetype =='extended_convex_hull':
                    toten=en*totatom
                else:
                    toten=en
                self.totalEnergy.append(toten)
            return self.totalEnergy
        else:
            for en,comp in zip(self.enthalpy_sub,self.compositions_vector_sub):
                totatom=sum(comp)
                #print(en)
                if self.filetype == 'extended_convex_hull':
                    toten = en * totatom
                else:
                    toten = en
                self.totalEnergy_sub.append(toten)
            return self.totalEnergy_sub


    def getHf(self): 
        # 如果文件形式是extended_convex_hull，输出以化学式为单位的整体能量；
        # 如果是individuals，则输出单个原子的平均能量
        if self.subsystem == None:
            for comp,toten in zip(self.compositions_vector,self.totalEnergy):
                Hf=toten
                #if self.filetype == 'extended_convex_hull':
                for elenum in range(len(comp)):
                    Hf-=comp[elenum]*self.atomenergy[elenum]
            '''else:
                    for elenum in range(len(comp)):
                        Hf-=comp[elenum]*self.atomenergy[elenum]
                    Hf = Hf/sum(comp)'''
                self.Hf.append(Hf)
        else:
            for comp,toten in zip(self.compositions_vector_sub,self.totalEnergy_sub):
                Hf = toten
                for elenum,atomen_sub in zip(range(len(comp)),self.subsystem):
                    Hf -= comp[elenum] * self.atomenergy[atomen_sub]
                self.Hf.append(Hf)
        return self.Hf #(ev/mole)
    
    
        

    def getRefs(self):
        # 获得ase画相图所需要的数据格式refs
        if self.subsystem == None:
            for comp,Hf in zip(self.compositions,self.Hf):
                ref=(comp,Hf)
                self.refs.append(ref)
        else:
            for comp,Hf in zip(self.compositions_sub,self.Hf):
                ref=(comp,Hf)
                self.refs.append(ref)
        return self.refs

    def remove_dup(self): 
        # 目前不对子系统进行筛选，只对全系统进行筛选
        # 内部函数，调用方便；提供了一个同名的外部函数，可以直接调用
        # 目前只针对id，compositions，compositions_vector,totalEnergy,hf,enthalpy 做筛选
        import numpy as np

        # dict={'id'=,'comp'=,'Hf'=}
        unduplicatedlist=[]
        unduplicatednum=[]
        num=0
        print('This might cost time, relax and take a cup of tea')
        for id,comp,Hf in zip(self.id,self.compositions_vector,self.Hf):
            totatom=sum(comp)
            Hf = Hf/totatom
            comp = np.array(comp)
            comp_normalized = comp / np.linalg.norm(comp)
            comp_dict = {'id': id, 'comp': comp,'Hf': Hf}
            if num == 0:
                unduplicatedlist.append(comp_dict)
            jud = 1
            for good_dict in unduplicatedlist:
                good_comp_normalized = good_dict['comp'] / np.linalg.norm(good_dict['comp'])
                if 1 - abs(np.dot(comp_normalized,good_comp_normalized)) < 1e-6 and Hf >= good_dict['Hf']:
                    jud -= 1
                    break
                elif 1 - abs(np.dot(comp_normalized,good_comp_normalized)) < 1e-6 and Hf < good_dict['Hf']:
                    jud -= 1
                    good_dict['id'] = id
                    good_dict['comp'] = comp
                    good_dict['Hf'] = Hf
                    break
                else:
                    jud -= 0
            if jud == 1:
                unduplicatedlist.append(comp_dict)
            num+=1
        for dic in unduplicatedlist:
            unduplicatednum.append(self.id.tolist().index(dic['id']))
        #print(unduplicatedlist[:5])
        remove_loc=[]
        for i in range(len(self.id)):
            if i not in unduplicatednum:
                remove_loc.append(i)
        self.id=np.delete(self.id,remove_loc)
        self.compositions=np.delete(self.compositions,remove_loc)
        self.compositions_vector=np.delete(self.compositions_vector,remove_loc,axis=0)
        self.enthalpy=np.delete(self.enthalpy,remove_loc)
        self.Hf=np.delete(self.Hf,remove_loc)
        self.totalEnergy=np.delete(self.totalEnergy,remove_loc)
        #print(self.compositions_vector[:10])
        return self
    
    def poscar_slice(self): 
        # 先阅读整个的extended_convex_hull_POSCARS
        poscars_ini = []
        for line in open(self.pospath,"r"):
            poscars_ini.append(line)
        poscarStart = 'EA' # USPEX输出的每个单独的POSCAR都以EA为开头
        #poscar_flies = [] # 一个列表，每一个元素就是一个单独的POSCAR文件
        start_loc = [] # 每个POSCAR开始的那一行的位置，保存为一个列表
        for line_num in range(len(poscars_ini)):
            if poscars_ini[line_num][0:2] == poscarStart:
                start_loc.append(line_num)
        #print(line_num)
        '''
        # 检查一下start_loc的长度是否和idlist的长度一致
        if len(start_loc) != len(self.id):
            print("error!amount of poscarfiles unequals to amount of ids in convexhull")
        # 检查一下每个POSCAR的ID是不是和convexhull.id里头的一致
        for num,id in zip(start_loc,self.id):
            if poscars_ini[num][0:7].strip() != poscarStart+str(id):
                print("error!POSCAR's id cannot correspond to ids in extendedconvexhull,please check!")
        '''
        # 将每个POSCAR文件切出来，存储到poscar_flies里头去
        # 这个循环里没法把最后一个poscar存进去，所以后面要单独写一条
        for i in range(len(start_loc)-1):
            poscar = []
            start = start_loc[i]
            end = start_loc[i+1]
            poscar = poscars_ini[start:end]
            self.poscar_silced.append(poscar)
        self.poscar_silced.append(poscars_ini[start_loc[-1]:])
        #self.poscar_silced = poscar_flies
    
    def poscar_selectByFitness(self):
        # 注意，这里筛选poscar的fitness依据是最开始设定的fitreq
        # 由于extend_convex_hull_POSCARS中的POSCAR与extended_convex_hull中的一致
        # 所以到这一步时，我们只需要前x个POSCAR就好了，这个x的值是经过fitreq,remove_dup等步骤筛选后的convexhull.id的长度
        for poscar in self.poscar_silced[0:len(self.id)]:
            #print(poscar)
            self.poscar_silced_selected.append(poscar)

    def poscar_buildfile(self):
        # 开始给每个POSCAR创建单独的文件夹，以方便后续VASP计算和吉布斯能量计算
        import os
        if os.path.exists(self.pos_savepath):
            print('save path already exists')
        else:
            os.mkdir(self.pos_savepath)
        # file_paths = [] 暂时好像没用
        poscarpaths = []
        if self.poscar_silced_selected != []:
            poscars = self.poscar_silced_selected
        else:
            poscars = self.poscar_silced
        for id,comp,poscar in zip(self.id,self.compositions,poscars):
            filename = 'EA'+str(int(id))+'_'+comp
            filepath = os.path.join(self.pos_savepath, filename)
            if os.path.exists(filepath) == False:
                os.mkdir(filepath)
            # file_paths.append(filepath) 暂时似乎没用
            poscarpath = os.path.join(filepath, 'POSCAR')
            poscarpaths.append(poscarpath)
            with open(poscarpath,'w') as f:
                for line in poscar:
                    f.write(line)
        # 保存POSCAR的路径，方便后续计算吉布斯能量和其他操作
        self.poscarpaths = poscarpaths
        
    def getPOTCAR(self):
        """
        参数解释
        potpath: POTCAR文件所在的路径
        pos_savepath: POSCAR文件所要存储的目标路径
        """
        import os
        from shutil import copyfile
        potcar_dic={} # 一共有多少种POTCAR，key为POTCAR文件名，value为POTCAR所含元素
        for potcarfile in os.listdir(self.potpath):
            if potcarfile[0:6] == 'POTCAR':
                potcar_dic[potcarfile]=potcarfile[7:]
        targetpaths = [] # 最终要把POTCAR都塞进这些文件夹里头
        for id,comp in zip(self.id,self.compositions):
            filename = 'EA'+str(int(id))+'_'+comp
            targetpath = os.path.join(self.pos_savepath, filename)
            #os.mkdir(filepath)
            targetpaths.append(targetpath)
        '''pot_ele = [] # 根据不同结构所含元素，指定所需要的POTCAR的种类'''
        for comp_vector,targetpath in zip(self.compositions_vector,targetpaths):
            ele = ''
            for i,j in zip(comp_vector,self.elements):
                if i != 0:
                    ele += j
            for key,value in potcar_dic.items():
                if ele == value:
                    potcarSource = os.path.join(self.potpath, key)
                    potcarTarget = os.path.join(targetpath, 'POTCAR') 
                    #print(potcarSource)
                    #print(potcarTarget)
                    copyfile(potcarSource,potcarTarget)
                    
    def getGibbs(self,mass,gels,T=[298.15]):
        # 参数解释
        # mass:元素质量，字典形式，key为元素符号，value为质量
        # gels:各元素的固体相吉布斯能，字典形式，key为元素符号，value为固体相吉布斯能
        # T:计算吉布斯能的温度，默认为298.15K，可自行设定，通常应该是一个列表
        #导入predcitG.py中的函数，计算吉布斯能,默认是计算298.15K的吉布斯能
        from PredictG import predictG

        #针对每个POSCAR文件，计算不同给定温度下的吉布斯能
        for poscarpath,compositon,hf in zip(self.poscarpaths, self.compositions,self.Hf):
            gibbs = predictG(
                compositon,
                hf,
                poscarpath,
                mass,
                gels,
            )
            dGs = []
            for t in T:
                dG = gibbs.dG(t, vol_per_atom=False)
                dGs.append(dG)
            # 以字典形式，存储comp和吉布斯能, key为comp，value为吉布斯能,作为列表元素，
            # 存入self.gibbs中
            gibbs_dict = {compositon:dGs}
            self.gibbs.append(gibbs_dict)


    
         
def remove_dup(ids,compositions,compositions_vectors,enthalpys,Hfs):
    import numpy as np

    # dict={'id'=,'comp'=,'Hf'=}
    unduplicatedlist=[]
    unduplicatednum=[]
    num=0
    for id,comp,Hf in zip(ids,compositions_vectors,Hfs):
        totatom=sum(comp)
        Hf = Hf/totatom
        comp = np.array(comp)
        comp_normalized = comp / np.linalg.norm(comp)
        comp_dict = {'id': id, 'comp': comp,'Hf': Hf}
        if num == 0:
            unduplicatedlist.append(comp_dict)
        jud = 1
        for good_dict in unduplicatedlist:
            good_comp_normalized = good_dict['comp'] / np.linalg.norm(good_dict['comp'])
            if 1 - abs(np.dot(comp_normalized,good_comp_normalized)) < 1e-6 and Hf >= good_dict['Hf']:
                jud -= 1
                break
            elif 1 - abs(np.dot(comp_normalized,good_comp_normalized)) < 1e-6 and Hf < good_dict['Hf']:
                jud -= 1
                good_dict['id'] = id
                good_dict['comp'] = comp
                good_dict['Hf'] = Hf
                break
            else:
                jud -= 0
        if jud == 1:
            unduplicatedlist.append(comp_dict)
        num+=1
    for dic in unduplicatedlist:
        unduplicatednum.append(ids.tolist().index(dic['id']))
    #print(unduplicatedlist[:5])
    remove_loc=[]
    for i in range(len(ids)):
        if i not in unduplicatednum:
            remove_loc.append(i)
    ids=np.delete(ids,remove_loc)
    compositions=np.delete(compositions,remove_loc)
    compositions_vectors=np.delete(compositions_vectors,remove_loc,axis=0)
    enthalpys=np.delete(enthalpys,remove_loc)
    Hfs=np.delete(Hfs,remove_loc)
    #print(self.compositions_vector[:10])
    return ids,compositions,compositions_vectors,enthalpys,Hfs

# 定义一个外部的函数，用于计算吉布斯能
def getGibbs(compositions,poscarpaths,Hfs,mass,gels,T=[298.15]):
    # 参数解释
    # mass:元素质量，字典形式，key为元素符号，value为质量
    # gels:各元素的固体相吉布斯能，字典形式，key为元素符号，value为固体相吉布斯能
    # T:计算吉布斯能的温度，默认为298.15K，可自行设定，通常应该是一个列表
    #导入predcitG.py中的函数，计算吉布斯能,默认是计算298.15K的吉布斯能
    from PredictG import predictG

    #针对每个POSCAR文件，计算不同给定温度下的吉布斯能
    gibbs = []
    for poscarpath,compositon,hf in zip(poscarpaths, compositions,Hfs):
        dGs = []
        for t in T:
            dG = predictG(
                compositon,
                hf,
                poscarpath,
                mass,
                gels,
            ).dG(t, vol_per_atom=False)
            dGs.append(dG)
        # 以字典形式，存储comp和吉布斯能, key为comp，value为吉布斯能,整体作为列表元素，添加进gibbs中
        gibbs_dict = {compositon:dGs}
        gibbs.append(gibbs_dict)
    return gibbs
