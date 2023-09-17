import sys

with open(sys.argv[1],'r') as f:
    with open(sys.argv[1]+'2','w') as outFile:
        for line in f:
            sp = line.strip().split()
            if len(sp) < 15:
                print('found a weird line in the RM output...',sp)
                continue

            elif '*' in line and len(sp) < 16:
                
                print('found a weird line in the RM output...',sp)
                continue

            else:
                outFile.write(line)

