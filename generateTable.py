'''
Generates a table containing all stentor genes, whether each gene belongs to an ortholog group, and the phyletic distributions of the group the gene is in.

Common abbreviations in script:
og : ortholog group
cg : ciliate-specific co-ortholog group
ng : 'NO_GROUP' - genes that matched individual proteins in OrthoMCL db that aren't members of groups
cam : ciliate/alveolate/metazoan
'''

from Bio import SeqIO
import csv

# dictionary of files containing ciliate/alveolate/metazoan ortholog group lists
cam_dict = {'tetrahymena_only.txt' : 'Ciliate-specific', 
            'api_only.txt' : 'Apicomplexan-specific',
            'alveolate_only.txt' : 'Alveolate-specific', 
            'CAM.txt' : 'Alveolates + Metazoa',
            'metazoa.txt' : 'Metazoa', 
            'metazoa_ciliates.txt' : 'Ciliates + Metazoa',
            'metazoa_apicomplexans.txt' : 'Apicomplexans + Metazoa'}

# dictionary of files containing domain ortholog group lists
domain_dict = {'bacteria_only.txt' : 'Bacteria Only', 
               'archaea_only.txt' : 'Archaea Only', 
               'euk_only.txt' : 'Eukaryotes Only', 
               'prokaryote_only.txt' : 'Bacteria + Archaea', 
               'bact_euk.txt' : 'Bacteria + Eukaryotes', 
               'arch_euk.txt' : 'Archaea + Eukaryotes', 
               'alldomains.txt' : 'All Domains'}

# OrthologGroup class definition
class OrthologGroup:
    def __init__(self, name, gene=''):
        self.name = name
        self.gene = gene
        self.cil = []
        self.cam = ''
        self.domain = ''

    def get_cil(self, t_bool, p_bool, o_bool):
        if t_bool:
            self.cil.append('Tetrahymena')
        if p_bool:
            self.cil.append('Paramecium')
        if o_bool:
            self.cil.append('Oxytricha')
        return

    def get_cam(self, cam_dict):
        if self.name.startswith('CiliateGroup') or self.gene.startswith('tthe|'):
            if self.cil:
                self.cam = 'Ciliate-specific'
            else:
                self.cam = 'None'
        elif self.name.startswith('Other'):
            if self.gene[:5] in ['mmus|', 'ecab|', 'mdom|', 'cbri|', 'tadh|', 'amel|', 'oana|', 'bmor|', 'cint|', 'apis|', 'cpip|',
                       'dmel|', 'rnor|', 'sman|', 'clup|', 'trub|', 'mmul|', 'isca|', 'ggal|', 'drer|', 'tnig|', 'nvec|']:
                if self.cil:
                    self.cam = 'Ciliates + Metazoa'
                else:
                    self.cam = 'Metazoa'
            elif self.gene[:5] in ['bbov|', 'ncan|', 'tgon|']:
                if self.cil:
                    self.cam = 'Alveolate-specific'
                else:
                    self.cam = 'Apicomplexan-specific'
            else:
                if self.cil:
                    self.cam = 'Ciliate-specific'
                else:
                    self.cam = 'None'
        else:
            for lst in cam_lists:
                if self.name in lst:
                    self.cam = lst[0]
                    break
        return
    
    def get_domain(self, domain_dict):
        if self.name.startswith('CiliateGroup'):
            self.domain = 'Eukaryotes Only'
        elif self.gene[:5] == 'cpne|':
            self.domain = 'Bacteria and Eukaryotes'
        elif self.name.startswith('Other'):
            self.domain = 'Eukaryotes Only'
        else:
            for lst in domain_lists:
                if self.name in lst:
                    self.domain = lst[0]
        return


# get list of Stentor genes
print 'Getting list of Stentor genes'
with open('../../AllGeneModels_Translated_2015_4_16.fasta', 'r') as genes_file:
    genes = [gene.id for gene in SeqIO.parse(genes_file, "fasta")]

gene_dict = {} #initializes dictionary mapping stentor genes to their ortholog groups
group_list = [] #initializes list of all ortholog groups


# assign OrthoMCL groups to genes
print 'Assigning OrthoMCL groups to genes'
with open('../stentor_orthomcl_results/orthologGroups', 'rb') as s_file:
    gene_groups = [[row[0], row[1], row[2]] for row in csv.reader(s_file, delimiter = '\t')]
for gene_group in gene_groups:
    if gene_group[1] != 'NO_GROUP': #curated OrthoMCL groups only
        if gene_group[1] not in [grp.name for grp in group_list]:
            group_list.append(OrthologGroup(gene_group[1]))
            gene_dict[gene_group[0]] = group_list[-1]
        else:
            gene_dict[gene_group[0]] = group_list[[grp.name for grp in group_list].index(gene_group[1])]
        

# assign Ciliate-specific co-ortholog groups to genes
print 'Assigning ciliate-specific co-ortholog groups to genes'

with open('../combined_orthomcl_results/paralogGroups', 'rb') as ciliate_file:
    ciliate_groups = list(csv.reader(ciliate_file, delimiter = '\t'))
for i in xrange(len(ciliate_groups)):
    ciliate_groups[i] = ['CiliateGroup' + str(i + 1)] + ciliate_groups[i]

with open('ciliatecoorthologs.txt', 'w') as ciliate_out: # print all ciliate groups to file
    ciliate_out.write('\n'.join(['\t'.join(group) for group in ciliate_groups]))

ciliate_groups_abbr = [[group[0]] + [item[0] for item in group[1:]] for group in ciliate_groups]
ciliate_dict = {} #initialize dictionary for assessing which species are in each ciliate ortholog group
for item in ciliate_groups_abbr:
    if 'g' in item:
        ciliate_dict[item[0]] = item[1:]

for ciliate_group in ciliate_groups:
    for item in ciliate_group:
        if item.startswith('g'):
            if ciliate_group[0] not in [grp.name for grp in group_list]:
                group_list.append(OrthologGroup(ciliate_group[0]))
                gene_dict[item] = group_list[-1]
            else:
                gene_dict[item] = group_list[[grp.name for grp in group_list].index(ciliate_group[0])]


# create groups for genes that have a "NO_GROUP" match
print 'Creating groups for genes that have a "NO_GROUP" match'

no_group = [row for row in gene_groups if row[1] == 'NO_GROUP']
no_group_genes = list(set([row[2] for row in no_group]))
other_groups = [['OtherGroup' + str(i+1), no_group_genes[i], []] for i in xrange(len(no_group_genes))]

for ng in no_group:
    for row in other_groups:
        if ng[2] == row[1]:
            ng.append(row[0])

for row in no_group:
    if row[-1] not in [grp.name for grp in group_list]:
        group_list.append(OrthologGroup(row[-1], row[2]))
        gene_dict[row[0]] = group_list[-1]
    else:
        gene_dict[row[0]] = group_list[[grp.name for grp in group_list].index(row[-1])]


# make note of genes that do not match any other gene
print 'Adding remaining genes with no matches'
for gene in genes:
    if gene not in gene_dict.keys():
        gene_dict[gene] = None
print '%d of %d genes added to dictionary' % (len(gene_dict.keys()), len(genes)) #checks that all stentor genes are in dictionary


# create lists of ortholog groups that are shared with Tetrahymena, Paramecium, or Oxytricha
with open('../tetrahymena_orthomcl_results/orthologGroups', 'rb') as t_file:
    t_ogs = list(set([row[1] for row in csv.reader(t_file, delimiter = '\t')]))
    t_nogroup = list(set([row[2] for row in csv.reader(t_file, delimiter = '\t') if row[1] == 'NO_GROUP']))
with open('../paramecium_orthomcl_results/orthologGroups', 'rb') as p_file:
    p_ogs = list(set([row[1] for row in csv.reader(p_file, delimiter = '\t')]))
    p_nogroup = list(set([row[2] for row in csv.reader(p_file, delimiter = '\t') if row[1] == 'NO_GROUP']))
with open('../oxytricha_orthomcl_results/orthologGroups', 'rb') as o_file:
    o_ogs = list(set([row[1] for row in csv.reader(o_file, delimiter = '\t')]))
    o_nogroup = list(set([row[2] for row in csv.reader(o_file, delimiter = '\t') if row[1] == 'NO_GROUP']))


# import lists of ortholog groups that are shared with ciliates, apicomplexans, metazoans, or a combination
cam_lists = []
for key in cam_dict.keys():
    with open(key, 'r') as cam_file:
        cam_lists.append([cam_dict[key]] + cam_file.read().split())


# import lists of ortholog groups that are shared with bacteria, archaea, eukaryotes, or a combination
domain_lists = []
for key in domain_dict.keys():
    with open(key, 'r') as domain_file:
        domain_lists.append([domain_dict[key]] + domain_file.read().split())


# get phyletic distribution info for each ortholog group, add as object attributes
print 'determining phyletic distributions of each ortholog group'
for grp in group_list:
    if grp.name.startswith('Ciliate'):
        grp.get_cil('T' in ciliate_dict[grp.name], 'P' in ciliate_dict[grp.name] or 'G' in ciliate_dict[grp.name], 'C' in ciliate_dict[grp.name])
    elif grp.name.startswith('Other'):
        grp.get_cil(grp.gene in t_nogroup, grp.gene in p_nogroup, grp.gene in o_nogroup)
    elif grp.name.startswith('OG5'):
        grp.get_cil(grp.name in t_ogs, grp.name in p_ogs, grp.name in o_ogs)
    else:
        print grp.name
    grp.get_cam(cam_dict)
    grp.get_domain(domain_dict)


# make table
print 'Writing table'
with open('stentor_orthology_table.txt', 'w') as outfile:
    for key in gene_dict.keys():
        if gene_dict[key]:
            outfile.write('%s\t%s\t%s\t%s\t%s\n' % (key, gene_dict[key].name, ', '.join(gene_dict[key].cil), gene_dict[key].cam, gene_dict[key].domain))
        else:
            outfile.write('%s\tNone\n' % key)
