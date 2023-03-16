def getSpaceGroup(path,symtole=1e-1):
    # path指存放POSCAR的文件夹目录的列表
    # 本脚本暂时无法处理子文件夹
    # symtole指找对称群时的tolerance
    import ase.io.vasp
    import ase.spacegroup
    import os
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
    print(sglist)
    return sglist

#print(getSpaceGroup(path=dirlist))