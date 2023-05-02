#!/bin/python
names=[] 
states=[]
unfinished_dic={} 
for line in open('./scan_check','r'):
    for i in range(len(line)):
        name=''
	state=''
	if line[i:i+5]=='state':
	    name=line[5:i]
	    state=line[i+6:]
	    names.append(name)
	    states.append(state)
for name,state in zip(names,states):
	if not ('reached' in state or 'writ' in state):
	    unfinished_dic[name]=state
lines=[]
for name,state in unfinished_dic.items():
    line='name='+name+'    '+'state='+state
    lines.append(line)
pattern=input('\npls choose a filter pattern,0-name and state,1-name:\n')
if pattern == 0:
    with open('./unfinished','w') as w:    
        for line in lines:
	    w.write(line)
else:
    with open('./name_unfinished','w') as w:
        for name in unfinished_dic.keys():
            name+='\n'
	    w.write(name)
