import csv

chip_file = 'data/Master_CHIP.csv'
pvalue_threshold = 1e-8

output_file = 'circos/chip_links.txt'
output_labels_file = 'circos/chip_labels.txt'


ch_genes = ['PRPF8', 'BRAF', 'JAK2', 'DNMT3A', 'SUZ12', 'PPM1D', 'ASXL1', 'GNB1', 'MPL', 'TET2', 'SRSF2', 'SF3B1', 'IDH2', 'CBL']

ensembl_gtf = '../chord_diagram/data/Homo_sapiens.GRCh38.107.gtf'
hgnc_file = '../chord_diagram/data/HGNC_names.txt'
olink_file = '../chord_diagram/data/olink_protein_map_1.5k_v1-Manifest.txt'

color = {'pos': 'positive', 'neg': 'negative', 'both': 'both'}

ensg2pos = {}
for line in open(ensembl_gtf):
	if line[:1] == '#':
		continue
	row = line.split()
	if row[2] == 'gene':
		start, end = int(row[3]), int(row[4])
		chrom = row[0]
		ensg = line.split('gene_id "')[1].split('"')[0]
		ensg2pos[ensg] = (chrom, start, end)

hgnc2ensg = {}
for x in csv.DictReader(open(hgnc_file), delimiter='\t'):
	hgnc2ensg[x['Approved symbol']] = x['Ensembl gene ID']
	for y in x['Previous symbols'].split(', '):
		if y not in hgnc2ensg:
			hgnc2ensg[y] = x['Ensembl gene ID']

prot2genes = {}
for row in csv.DictReader(open(olink_file), delimiter='\t'):
	acc = row['UKBPPP_ProteinID'].split(':')[0]
	prot2genes[acc] = prot2genes.get(acc, []) + [row['HGNC.symbol']]
	hgnc2ensg[row['HGNC.symbol']] = row['ensembl_id']

complexes = set()
interactions = {}
for row in csv.DictReader(open(chip_file), delimiter=';'):
	if 'VAF10' in row['geno']:
		gene = row['geno'].split('_')[0]

		if gene != 'TET2':
			continue

		if gene == 'CH':
			genes = ch_genes
		else:
			genes = [gene]



		proteins = list(set(prot2genes[row['pheno'].split(':')[0]]))

		for gene in genes:
			for prot in proteins:
				p = float(row['p'].replace(',','.'))
				beta = float(row['beta'].replace(',','.'))
				inter = 'neg' if beta < 0 else 'pos'

				if p <= pvalue_threshold:
					if (gene, prot) in interactions:
						if interactions[(gene, prot)] != inter:
							interactions[(gene, prot)] = 'both'
					else:
						interactions[(gene, prot)] = inter


print(len(interactions), 'total interactions')
print()

genes_all = set()
with open(output_file, 'w') as f:
	for gene, prot in interactions:
		gene_pos, prot_pos = ensg2pos[hgnc2ensg[gene]], ensg2pos[hgnc2ensg[prot]]
		genes_all.add(gene)
		genes_all.add(prot)
		
		source_chrom = 'hs{}'.format(gene_pos[0])
		source_start = gene_pos[1]
		source_end = gene_pos[2]#source_start + gene_width

		target_chrom = 'hs{}'.format(prot_pos[0])
		target_start = prot_pos[1]
		target_end = prot_pos[2]#target_start + gene_width

		print(source_chrom, source_start, source_end, target_chrom, target_start, target_end, 
			'color='+color[interactions[(gene, prot)]], 
			sep='\t', file=f)

with open(output_labels_file, 'w') as f:
	for gene in genes_all:

		gene_pos = ensg2pos[hgnc2ensg[gene]]
		
		source_chrom = 'hs{}'.format(gene_pos[0])
		source_start = gene_pos[1]
		source_end = gene_pos[2]

		print(source_chrom, source_start, source_end, gene, sep='\t', file=f)