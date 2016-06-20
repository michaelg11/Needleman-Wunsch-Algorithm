def zeros(m,n):
	return [[0]*n for _ in xrange(m)]
	
def match_score(k,l):
	fname='/Users/mridulgarg11/Documents/Spring\'16/ECEN_689/BLOSUM62.txt'
	handle=open(fname)
	matrix=[]
	col,row={},{}
	idx,idr=0,0
	first_line=0
	for line in handle:
		if not line.startswith('#'):
			if first_line==0:
				first_line=1
				for i in line.split():
					col[i]=idx
					idx+=1
			else:
				a=line.split()
				x=a.pop(0)
				row[x]=idr
				matrix.append(a)
				idr+=1
	return int(matrix[col[k]][row[l]])

def nw (seq1,seq2,pen):
	m=len(seq1)
	n=len(seq2)
	table=zeros(m+1,n+1)
	
	for i in range(0,m+1):
		table[i][0]=pen*i
	for j in range(0,n+1):
		table[0][j]=pen*j
	for i in range(1,m+1):
		for j in range(1,n+1):
			match=table[i-1][j-1]+match_score(seq1[i-1],seq2[j-1])
			delete=table[i-1][j]+pen
			insert=table[i][j-1]+pen
			table[i][j]=max(match,insert,delete)
		
		
	align1,align2='',''
	i,j=m,n
	while i>0 and j>0:
		current=table[i][j]
		diaganol=table[i-1][j-1]
		up=table[i][j-1]
		left=table[i-1][j]
		
		if current==diaganol+match_score(seq1[i-1],seq2[j-1]):
			align1+=seq1[i-1]
			align2+=seq2[j-1]
			i,j=i-1,j-1
			
		elif current==left+pen:
			align1+=seq1[i-1]
			align2+='-'
			i-=1
		
		elif current==up+pen:
			align1+='-'
			align2+=seq2[j-1]
			j-=1
		
	while i>0:
		align1+=seq1[i-1]
		align2+='-'
		i-=1			
	while j>0:
		align1+='-'
		align2+=seq2[j-1]
		j-=1
		
	seqalign(align1,align2)
	
def seqalign (align1,align2):
	align1=align1[::-1]
	align2=align2[::-1]
	i,j=0,0
	
	seqscore=0
	identity=0
	symbol=''
	
	for i in range(0,len(align1)):
		if align1[i]==align2[i]:
			identity+=1
			seqscore+=match_score(align1[i],align2[i])
			symbol+= '|'
			
		elif align1[i]=='-' or align2[i]=='-':
			seqscore+=pen
			symbol+=' '
		
		elif align1[i]!=align2[i] and align1[i]!='-' and align2[i]!='-':
			seqscore+=match_score(align1[i],align2[i])
			symbol+=' '
				
	identity= float(identity)/len(align1) *100
	
	print 'Sequence similarity:', seqscore
	print 'Identity Score', identity,'%'
	print align1
	print symbol
	print align2
	

seq4='ADCNSRQCLCRPM'
seq5='ASCSNRCKCRDP'
seq1='GLSDGEWQLVLNVWGKVEADVAGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGNTVLTALGGILKKKGHHEAELTPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMSKALELFRNDMAAKYKELG'
seq2='GLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELG'
seq3='GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDRFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELG'
pen=-8

nw(seq1,seq2,pen)
nw(seq1,seq3,pen)
