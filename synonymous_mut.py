
# <--------------------DEFENITIONS-------------------->

#Reading fasta sequence from file and storing sequence without header in a string
def read_fasta(file1):
	lines = []
	for a in open (file1):
		a = a.strip()
		if a != "":
			lines.append(a)
	
	lines.pop(0)	
	seq = ''.join(lines)
	return seq


#Reading a CSV file and creating dictionary of SNPs position and variant
def read_poly (file2):
	pairs = {}
	for b in open (file2):
		b = b.strip()
		if b != "":
			fields = b.split(',')
			pairs[fields[0]] = fields[11]	
	del pairs[list(pairs)[0]]
	return pairs


#translating genomic region using standard codon table
def comp_translation (wild,mut):
	codons = {'TTT':'F','TCT':'S','TAT':'Y','TGT':'C','TTC':'F','TCC':'S','TAC':'Y','TGC':'C','TTA':'L','TCA':'S','TAA':'*','TGA':'*','TTG':'L','TCG':'S','TAG':'*','TGG':'W',
	'CTT':'L','CCT':'P','CAT':'H','CGT':'R','CTC':'L','CCC':'P','CAC':'H','CGC':'R','CTA':'L','CCA':'P','CAA':'Q','CGA':'R','CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
	'ATT':'I','ACT':'T','AAT':'N','AGT':'S','ATC':'I','ACC':'T','AAC':'N','AGC':'S','ATA':'I','ACA':'T','AAA':'K','AGA':'R','ATG':'M','ACG':'T','AAG':'K','AGG':'R',
	'GTT':'V','GCT':'A','GAT':'D','GGT':'G','GTC':'V','GCC':'A','GAC':'D','GGC':'G','GTA':'V','GCA':'A','GAA':'E','GGA':'G','GTG':'V','GCG':'A','GAG':'E','GGG':'G'}
	wild_peptide=[]
	mut_peptide=[]
	#translating wild and mutant region
	for i in range(0,len(mut),3):
		wild_peptide.append(codons[''.join(wild[i:i+3])])
		mut_peptide.append(codons[''.join(mut[i:i+3])])
	
	#comparing peptides and determining position of difference
	count = 1
	for w,m in zip(wild_peptide,mut_peptide):
		if w != m:
			replacement = w+" -> "+m+" at "+str(count)+" ("+str(len(wild_peptide))+"aa long)"
			break
		else:
			replacement = "Synonymous"
		count = count + 1
	return replacement
	

#getting regions for translation with SNP information and reporting changes in peptide
def trans_region (start,end,pos,var):
	mut_genome = list(genome)
	mut_genome[pos] = var
	wild_region = list(genome [start:end])      #region from the wild type genome
	mut_region = mut_genome [start:end]         #region from the mutant genome
	change = comp_translation(wild_region,mut_region)
	return change 
	
	
# <----------------------PROGRAM---------------------->


#getting genome from fasta file
fasta_file_name = "ncbi_dengu1.fasta"
genome = read_fasta (fasta_file_name)

#getting dictionary of SNPs from CSV file
poly_file_name = "poly_dengu1.csv"
poly = read_poly (poly_file_name)

#mature peptide coding regions
#start,end,gene,locus_tag,gene_synonym,product,protein_id

peptide_regions = [[95,436,"POLY","DV1_gp1","polyprotein gene","anchored capsid protein ancC","NP_722457.2"],
[95,394,"POLY","DV1_gp1","polyprotein gene","capsid protein C","NP_722466.2"],
[437,934,"POLY","DV1_gp1","polyprotein gene","membrane glycoprotein precursor prM","NP_733807.2"],
[437,709,"POLY","DV1_gp1","polyprotein gene","protein pr","YP_009164956.1"],
[710,934,"POLY","DV1_gp1","polyprotein gene","membrane glycoprotein M","NP_722459.2"],
[935,2419,"POLY","DV1_gp1","polyprotein gene","envelope protein E","NP_722460.2"],
[2420,3475,"POLY","DV1_gp1","polyprotein gene","nonstructural protein NS1","NP_722461.1"],
[3476,4129,"POLY","DV1_gp1","polyprotein gene","nonstructural protein NS2A","NP_733808.1"],
[4130,4519,"POLY","DV1_gp1","polyprotein gene","nonstructural protein NS2B","NP_733809.1"],
[4520,6376,"POLY","DV1_gp1","polyprotein gene","nonstructural protein NS3","NP_722463.1"],
[6377,6757,"POLY","DV1_gp1","polyprotein gene","nonstructural protein NS4A","NP_733810.1"],
[6758,6826,"POLY","DV1_gp1","polyprotein gene","protein 2K","NP_722467.1"],
[6827,7573,"POLY","DV1_gp1","polyprotein gene","nonstructural protein NS4B","NP_733811.1"],
[7574,10270,"POLY","DV1_gp1","polyprotein gene","RNA-dependent RNA polymerase NS5","NP_722465.1"]]

		
#iterating over peptide region and cheking SNPs in those region

print ("protein_id","gene","product\t\t\t\t","pos","target","variant","mutation type",sep = "\t")
for reg in peptide_regions:
	flag = 0 #flag if no SNP found in region
	for k,v in poly.items():
		if int(k)>=(reg[0]-1) and int(k)<=(reg[1]-1):
			pep_match_result = trans_region(reg[0]-1,reg[1],int(k)-1,v)
			print (reg[6],reg[2],reg[5],k,genome[int(k)-1],v,pep_match_result,sep = "\t")
			flag = 1
	if flag == 0:
		print (reg[6],reg[2],reg[5],"\tNo SNP",sep ="\t")
	
	
	
