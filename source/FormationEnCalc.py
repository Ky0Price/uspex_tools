def FormationEnCalc(idlist,system,moleEn,atomEn,op=1):
    # idlist是每一个结构对应的一个ID，后续会用到，随便什么形式，但在本脚本中应统一
    # moleEn指分子总能量，是一个列表
    # system指体系的组成向量,形式为'Ba3_C1_N0'或'Ba3 C1 N0',总是每个元素和前一个元素+个数的组合之间需要有间隔符
    # atomEn指的是组成分子前每个原子的基态能量，需要自己手动写个字典传入，比如{'Li':-2,'C':-9}
    # 形成焓计算公式deltaE(mAnB) = E(mAnB)-m*E(A)-n*E(B)
    # op选择处理模式，1会输出分子的总生成焓，适用于ase作相图，0会输出分子中每个原子的平均生成焓
    FElist=[] #存放最终的ID-组分-生成焓的对应字典的列表
    for id,comp,en in zip(idlist,system,moleEn):
        # 计算体系中原子种类数和各原子种类对应的原子数目
        FEdict = {}
        FormationEn = en
        complist = comp.split('_') # 按元素以你指定的间隔符切割一下，所以要求组分形式为'Ba3_C1_N0'此类，也可以用空格'Ba3 C1 N0'，则split中的参数要改为(' ')。
        totalAtomNumber = 0
        composition = ''  # 用于生成方便输入ASE的组成格式
        for ele in complist:
            element = ''  # 记录元素种类
            number = ''  # 记录原子个数
            for i in ele:
                if i.isalpha():
                    element += i
                else:
                    number += i
            totalAtomNumber += eval(number)
            composition += element+number
            FormationEn -= atomEn[element]*eval(number) # 减去对应原子的基态能量
        #print(composition)
        if op==1:
            FEdict['composition']=composition 
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
    



