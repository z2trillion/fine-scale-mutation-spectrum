group_to_populations = {
    'EAS': ['CHB','JPT','CHS','CDX','KHV','CHD'],
    'EUR': ['CEU','TSI','GBR','FIN','IBS'],
    'AFR': ['YRI','LWK','GWD','MSL','ESN', 'ACB', 'ASW'],
    'SAS': ['GIH','PJL','BEB','STU','ITU'],
    'AMR': ['CLM','MXL','PUR','PEL'],
}

groups = group_to_populations.keys()

population_to_group = {}
for group, populations in group_to_populations.iteritems():
    for population in populations:
        population_to_group[population]=group

populations = population_to_group.keys()

sample_id_to_population = {}
with open('data/1000genomes_phase3_sample_IDs.txt') as sample_id_lines:
    for line in sample_id_lines:
        sample_id, _, population, _ = line.split(None, 3)
        sample_id_to_population[sample_id] = population
