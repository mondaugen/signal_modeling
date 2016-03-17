import sys
ptable=dict()
# to match MIDI pitch lowest note is 12
n=12
for m in xrange(0,11):
    for note in ['c','d','e','f','g','a','b']:
        ptable[note+str(m)]=n
        n=n+1
        if note != 'e' and note != 'b':
            ptable[note+'#'+str(m)]=n
            n=n+1

for line in sys.stdin:
    line=line.strip().lower()
    print line+' '+str(ptable[line])
