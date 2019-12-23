codons = {'TTT':'F','TCT':'S','TAT':'Y','TGT':'C','TTC':'F','TCC':'S','TAC':'Y','TGC':'C','TTA':'L','TCA':'S','TAA':'*','TGA':'*','TTG':'L','TCG':'S','TAG':'*','TGG':'W',
'CTT':'L','CCT':'P','CAT':'H','CGT':'R','CTC':'L','CCC':'P','CAC':'H','CGC':'R','CTA':'L','CCA':'P','CAA':'Q','CGA':'R','CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
'ATT':'I','ACT':'T','AAT':'N','AGT':'S','ATC':'I','ACC':'T','AAC':'N','AGC':'S','ATA':'I','ACA':'T','AAA':'K','AGA':'R','ATG':'M','ACG':'T','AAG':'K','AGG':'R',
'GTT':'V','GCT':'A','GAT':'D','GGT':'G','GTC':'V','GCC':'A','GAC':'D','GGC':'G','GTA':'V','GCA':'A','GAA':'E','GGA':'G','GTG':'V','GCG':'A','GAG':'E','GGG':'G'}

records = []
for n in open (file2):
	n = n.strip().replace('-','0')
	records.append(n)
records.pop(0)

pairs = {}
for b in records:
	if b != "":
		fields = b.split(',')
		if fields[1] != '0':
			variant = 'A'
			var_depth = int(fields[7])
		elif fields[2] != '0':
			variant = 'C'
			var_depth = int(fields[8])
		elif fields[3] != '0':
			variant = 'G'
			var_depth = int(fields[9])
		elif fields[4] != '0':
			variant = 'T'
			var_depth = int(fields[10])
		
		if fields[5] == 'a':
			target = 'A'
			tar_depth = int(fields[7])
		elif fields[5] == 'c':
			target = 'C'
			tar_depth = int(fields[8])
		elif fields[5] == 'g':
			target = 'G'
			tar_depth = int(fields[9])
		elif fields[5] == 't':
			target = 'T'
			tar_depth = int(fields[10])
	pairs[fields[0]] = [variant,var_depth,target,tar_depth]
return pairs
