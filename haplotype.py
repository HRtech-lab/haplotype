import pandas as pd
import numpy as np
import pysam
import allel
import re
import networkx as nx
import copy

# function to remove characters from string
def remove_character(s):
  return re.sub(r"[^0-9]","",s)

# function to find first character index in string
def find_first_letter_index(st):
  res = None
  temp = re.search(r'[a-z]', st, re.I)
  if temp is not None:
      res = temp.start()
  return res

# merging duplicate patterns data
def better_np_unique(arr):
    sort_indexes = np.argsort(arr)
    arr = np.asarray(arr)[sort_indexes]
    vals, first_indexes, inverse, counts = np.unique(arr,
        return_index=True, return_inverse=True, return_counts=True)
    indexes = np.split(sort_indexes, first_indexes[1:])
    for x in indexes:
        x.sort()
    return vals, indexes, inverse, counts

# function to add index into index corresponding value list
def add_index(li,ind):
  temp=[]
  for l,i in zip(li,ind):
    temp.append(l+i)
  return temp

# getting merged heplotype 1 list including corresponding indexes
def merge_duplicates(li):
  g = nx.Graph()
  ipath = li
  temp=[]
  for p in ipath:
    g.add_edges_from(zip(p, p[1:]))
  for c in nx.connected_components(g):
    temp.append(list(c))
  return temp

def get_unique_value_index_pair(li, index):
  li_arr = np.array(li)
  target_list = li_arr[:,index].tolist()
  uniqueChrom, indexes, _, _ = better_np_unique(target_list)
  paired_data = zip(uniqueChrom, indexes)
  return paired_data

def get_indexed_dictionary(li, li_paired_data):
  dpm_dict ={}
  for val in li_paired_data:
    temp = []
    for i in val[1]:
      temp.append(li[i])
    dpm_dict[val[0]] = temp
  return dpm_dict

#Working directory path
base_path =''
#gencode = pd.read_table(base_path + "alphonso.gtf", comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])
#Now extracting genes information from our dataframe
#gencode_genes_list = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1).values.tolist()

#Function to convert vcf entries into dataframe
import io
import os
import pandas as pd


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
  
 
vcf_df = read_vcf("napus.vcf")
#Getting information of the sample column
vcf_df = vcf_df.set_axis([*vcf_df.columns[:-1], 'add'], axis=1, inplace=False)
#Filtering heterozygous information from the vcf file
vcf_df = vcf_df[vcf_df['add'].str.contains('0/1', regex=False, case=False, na=False)]
vcf_dict = {}
# creating chrom indexed vcf dictionary in vcf_dict
uniqueChrom = vcf_df.CHROM.unique()
for entry in uniqueChrom:
  data = vcf_df.loc[(vcf_df['CHROM']  == entry)]
  vcf_dict[entry]={}
  for index, row in data.iterrows():
    if len(row['ALT']) == 1 and len(row['REF']) == 1:
      vcf_dict[entry][row['POS']] = [row['ALT'], row['REF'],row["add"]]
      
      
# samfile reader
samfile = pysam.AlignmentFile(base_path + "napus_sorted.bam","rb")
# creating dictionary of vcf compared patterns data
vcf_compared_data = []
# vcd = open(base_path + "wed_00BEFORE_vcf_compared_data.txt", "a")
# gathering read mate patterns in samfile_dict
j = 0
for g in vcf_dict.keys():
  vcf_chrom_dict = vcf_dict[g]
  bam_list = [*samfile.fetch(g)]
  bam_list.sort(key=lambda x: x.qname)
  i = 0
  while i < (len(bam_list)-1):
    rm_list = []
    if bam_list[i].is_proper_pair == True and bam_list[i].is_secondary == False and bam_list[i].mapq >= 20:
      remove_from = 0
      seq_start = 0
      i_value = 0
      if bam_list[i].cigarstring:  #Dealing NONE cases
        fli = find_first_letter_index(bam_list[i].cigarstring)
        # Handling cases if S or I in the first letter
        if  bam_list[i].cigarstring[fli] == "S":
          seq_start = int(bam_list[i].cigarstring[0:fli])
        elif  bam_list[i].cigarstring[fli] == "I":
          i_value = int(bam_list[i].cigarstring[0:fli])

        # Handling cases if I present but not first as character
        i_index = bam_list[i].cigarstring[:-1].find("I")
        if i_index != -1 and bam_list[i].cigarstring[fli] != "I":
          if  bam_list[i].cigarstring[fli] == "M":
            remove_from = int(bam_list[i].cigarstring[0:fli])
          loop_counter = i_index
          i_number = ""
          while loop_counter > fli:
            if bam_list[i].cigarstring[loop_counter-1].isdigit():
              i_number = i_number + bam_list[i].cigarstring[loop_counter-1]
              loop_counter = loop_counter - 1
            else:
              break
          i_value = int(i_number[::-1])

        for point in zip(bam_list[i].positions, bam_list[i].seq[seq_start:][:remove_from] + bam_list[i].seq[seq_start:][remove_from+i_value:]):
          if (point[0]+1 in vcf_chrom_dict.keys()) and (point[1] in vcf_chrom_dict[point[0]+1]):
            rm_list.append(str(point[0]+1)+point[1])

    if bam_list[i].qname == bam_list[i+1].qname and bam_list[i+1].is_proper_pair == True and bam_list[i+1].is_secondary == False and bam_list[i+1].mapq >= 20:
      remove_from = 0
      seq_start = 0
      i_value = 0
      if bam_list[i+1].cigarstring:
        fli = find_first_letter_index(bam_list[i+1].cigarstring)
        # Handling cases if S or I in the first letter
        if  bam_list[i+1].cigarstring[fli] == "S":
          seq_start = int(bam_list[i+1].cigarstring[0:fli])
        elif  bam_list[i+1].cigarstring[fli] == "I":
          i_value = int(bam_list[i+1].cigarstring[0:fli])

        # Handling cases if I present but not first as character
        i_index = bam_list[i+1].cigarstring[:-1].find("I")
        if i_index != -1 and bam_list[i+1].cigarstring[fli] != "I":
          if  bam_list[i+1].cigarstring[fli] == "M":
            remove_from = int(bam_list[i+1].cigarstring[0:fli])
          loop_counter = i_index
          i_number = ""
          while loop_counter > fli:
            if bam_list[i+1].cigarstring[loop_counter-1].isdigit():
              i_number = i_number + bam_list[i+1].cigarstring[loop_counter-1]
              loop_counter = loop_counter - 1
            else:
              break
          i_value = int(i_number[::-1])

        for point in zip(bam_list[i+1].positions, bam_list[i+1].seq[seq_start:][:remove_from] + bam_list[i+1].seq[seq_start:][remove_from+i_value:]):
          if (point[0]+1 in vcf_chrom_dict.keys()) and (point[1] in vcf_chrom_dict[point[0]+1]):
            rm_list.append(str(point[0]+1)+point[1])
        # Incrementing i if read and mate both exist
        i = i+1

    # removing duplicates from rm_list
    rm_list = list(dict.fromkeys(rm_list))
    if len(rm_list) > 0:
      rm_wc_list = list(map(remove_character, rm_list))
      if len(rm_wc_list) == len(list(dict.fromkeys(rm_wc_list))):
        vcf_compared_data.append([bam_list[i].qname, g, rm_list])
        # vcd.write(bam_list[i].qname+";"+g+";"+str(rm_list)+"\n")

    # loop increment
    i = i+1

  # Loop breaking
  j = j+1
  print(j, end=" ")
#   if j > 0:
#     break
    
# vcd.close()


#creating chrom indexed dictionary
unique_chrom_pairs = get_unique_value_index_pair(vcf_compared_data, 1)
chrom_indexed_dictionary = get_indexed_dictionary(vcf_compared_data, unique_chrom_pairs)

# Merging chrom wise dito duplicate patterns according to count 
for chrom in chrom_indexed_dictionary.keys():
  chrom_data = chrom_indexed_dictionary[chrom]
  vals, indexes, _, counts = better_np_unique(np.array(chrom_data, dtype=object)[:,2])
  temp = zip(vals, indexes,counts)
  duplicate_merged = []
  for val in temp:
    duplicate_merged.append([val[0], val[2]])
  chrom_indexed_dictionary[chrom] = duplicate_merged
  
  fmo = open(base_path + "Test_first_merging_outputs.txt", "a")

#print after sorting w.r.t chrom and genes
for chrom in sorted([*chrom_indexed_dictionary.keys()], key=lambda x: int(remove_character(x))):
  for entry in chrom_indexed_dictionary[chrom]:
    fmo.write( "["+str(chrom)+";"+str(entry[0])+";"+str(entry[1])+"]\n" )

fmo.close()

# Removing duplicate positions with same characters if there is any same position with different alt from pattern having less count (chrom wise)
for chrom in chrom_indexed_dictionary.keys():
  data = chrom_indexed_dictionary[chrom]
  rm_list = copy.deepcopy(np.array(data, dtype=object)[:,0].tolist())
  rm_merged_list = merge_duplicates( add_index(rm_list, [[item] for item in range(len(rm_list))]) )
  for i in range(len(rm_merged_list)):
    items = rm_merged_list[i]
    indexes = [e for e in items if isinstance(e, int)]
    j = 0
    while j < len(indexes):
      rm_pairs = rm_list[indexes[j]]
      rm_pairs_without_alt = list(map(remove_character, rm_pairs))
      rm_count = chrom_indexed_dictionary[chrom][indexes[j]][1]
      for k in range(j+1, len(indexes)):
        rm_pairs_2 = rm_list[indexes[k]]
        rm_pairs_without_alt_2 = list(map(remove_character, rm_pairs_2))
        rm_count_2 = chrom_indexed_dictionary[chrom][indexes[k]][1]
        cwa = set(rm_pairs).intersection(rm_pairs_2)
        cwoa = set(rm_pairs_without_alt).intersection(rm_pairs_without_alt_2)
        if len(cwoa) > len(cwa) and len(cwa) > 0:
          if rm_count < rm_count_2:
            for entry in cwa:
              rm_pairs.remove(entry)
          else:
            for entry in cwa:
              rm_pairs_2.remove(entry)
      chrom_indexed_dictionary[chrom][indexes[j]][0] = copy.deepcopy(rm_pairs)
      j = j+1
      
   # merging the hyplotypes having atleast 1 same value
for chrom in chrom_indexed_dictionary.keys():
  chrom_data = chrom_indexed_dictionary[chrom]
  rm_list = np.array(chrom_data, dtype=object)[:,0].tolist()
  rm_merged_list = merge_duplicates( add_index(rm_list, [[item] for item in range(len(rm_list))]) )
  merged_patterns = []
  for i in range(len(rm_merged_list)):
    items = rm_merged_list[i]
    indexes = [e for e in items if isinstance(e, int)]
    rm_pairs = []
    for index in indexes:
      cwa = len(set(chrom_data[index][0]).intersection(rm_pairs))
      cwoa = len(set(list(map(remove_character, chrom_data[index][0]))).intersection(list(map(remove_character, rm_pairs))))
      if cwa == cwoa:
        rm_pairs = rm_pairs + chrom_data[index][0]
    merged_patterns.append(sorted(list(dict.fromkeys(rm_pairs)), key=lambda x: int(remove_character(x))))
    
    
    
    
    # removing ALT form all enteries of rm patterns so that we can find haplotypes of each other
alt_removed_dict = {}
for chrom in chrom_indexed_dictionary:
  alt_removed_dict[chrom] = {}
  chrom_data = chrom_indexed_dictionary[chrom]
  collect = []
  for entry in chrom_data:
      collect.append(list(map(remove_character, entry)))
  alt_removed_dict[chrom] = collect

# finding haplotypes in which first index contain the largest pattern
for chrom in alt_removed_dict.keys():
  rm_wc_list = alt_removed_dict[chrom]
  rm_merged_list = merge_duplicates( add_index(rm_wc_list, [[item] for item in range(len(rm_wc_list))]) )
  haplotype_merged = []
  for i in range(len(rm_merged_list)):
    items = rm_merged_list[i]
    indexes = sorted([e for e in items if isinstance(e, int)], key=lambda x: len(rm_wc_list[x]), reverse = True)
    visited = []
    j = 0
    while j < len(indexes):
      if not (j in visited):
        visited.append(j)
        haplotypes = []
        haplotypes.append(chrom_indexed_dictionary[chrom][indexes[j]])
        large_pattern = rm_wc_list[indexes[j]]
        for k in range(j+1, len(indexes)):
          small_pattern = rm_wc_list[indexes[k]]
          if set(small_pattern) <= set(large_pattern):
            haplotypes.append(chrom_indexed_dictionary[chrom][indexes[k]])
            visited.append(k)
        haplotype_merged.append(copy.deepcopy(haplotypes))
      j = j+1
  chrom_indexed_dictionary[chrom] = copy.deepcopy(haplotype_merged)
  
  # creating chrom based haplotype dictionary
chrom_based_dict_of_haplotypes = {}
for key in chrom_indexed_dictionary.keys():
  data = chrom_indexed_dictionary[key]
  haplotypes_list = []
  for hyplotype_list in data:
    first_haplotype = hyplotype_list[0]
    haplotypes_list.append(first_haplotype)
    second_haplotype = []
    for k in range(len(first_haplotype)):
      position = int(remove_character(first_haplotype[k]))
      haplotype_1 = first_haplotype[k][-1]
      if vcf_dict[key][position][0] == haplotype_1:
        second_haplotype.append(str(position)+vcf_dict[key][position][1])
      elif vcf_dict[key][position][1] == haplotype_1:
        second_haplotype.append(str(position)+vcf_dict[key][position][0])
    haplotypes_list.append(second_haplotype)
  chrom_based_dict_of_haplotypes[key] = haplotypes_list
  
  # Performing last merging to merge any remaining hyplotypes having atleast 1 same value
for chrom in chrom_based_dict_of_haplotypes.keys():
  chrom_data = chrom_based_dict_of_haplotypes[chrom]
  rm_merged_list = merge_duplicates( add_index(chrom_data, [[item] for item in range(len(chrom_data))]) )
  merged_patterns = []
  for i in range(len(rm_merged_list)):
    items = rm_merged_list[i]
    indexes = [e for e in items if isinstance(e, int)]
    rm_pairs = []
    for index in indexes:
      cwa = len(set(chrom_data[index]).intersection(rm_pairs))
      cwoa = len(set(list(map(remove_character, chrom_data[index]))).intersection(list(map(remove_character, rm_pairs))))
      if cwa == cwoa:
        rm_pairs = rm_pairs + chrom_data[index]
    merged_patterns.append(sorted(list(dict.fromkeys(rm_pairs)), key=lambda x: int(remove_character(x))))
  chrom_based_dict_of_haplotypes[chrom] = merged_patterns

  # removing ALT form all enteries of rm patterns so that we can find haplotypes of each other
alt_removed_dict = {}
for chrom in chrom_based_dict_of_haplotypes:
  alt_removed_dict[chrom] = {}
  chrom_data = chrom_based_dict_of_haplotypes[chrom]
  collect = []
  for entry in chrom_data:
      collect.append(list(map(remove_character, entry)))
  alt_removed_dict[chrom] = collect

# finding haplotypes in which first index contain the largest pattern
for chrom in alt_removed_dict.keys():
  rm_wc_list = alt_removed_dict[chrom]
  rm_merged_list = merge_duplicates( add_index(rm_wc_list, [[item] for item in range(len(rm_wc_list))]) )
  haplotype_merged = []
  for i in range(len(rm_merged_list)):
    items = rm_merged_list[i]
    indexes = sorted([e for e in items if isinstance(e, int)], key=lambda x: len(rm_wc_list[x]), reverse = True)
    visited = []
    j = 0
    while j < len(indexes):
      if not (j in visited):
        visited.append(j)
        haplotypes = []
        haplotypes.append(chrom_based_dict_of_haplotypes[chrom][indexes[j]])
        large_pattern = rm_wc_list[indexes[j]]
        for k in range(j+1, len(indexes)):
          small_pattern = rm_wc_list[indexes[k]]
          if set(small_pattern) <= set(large_pattern):
            haplotypes.append(chrom_based_dict_of_haplotypes[chrom][indexes[k]])
            visited.append(k)
        haplotype_merged.append(copy.deepcopy(haplotypes))
      j = j+1
  chrom_based_dict_of_haplotypes[chrom] = copy.deepcopy(haplotype_merged)
  
  smo = open(base_path + "test_Second_merging_HASSAN.txt", "a")
#print after sorting w.r.t chrom and genes with tabed spaces
i = 1
j = 1
for key in sorted([*chrom_based_dict_of_haplotypes.keys()], key=lambda x: int(remove_character(x))):
  data = sorted(chrom_based_dict_of_haplotypes[key], key=lambda x: int(remove_character(x[0][0])))
  for hyplotypes_list in data:
    reverse_flag = False
    first_entry = hyplotypes_list[0]
    if len(first_entry)>1:
        span = int(remove_character(first_entry[len(first_entry)-1])) - int(remove_character(first_entry[0]))
        smo.write("BLOCK: offset: "+str(j)+" len: "+str(len(first_entry))+" phased: "+str(len(first_entry))+" SPAN: "+str(span)+"\n")
        pos_1 = int(remove_character(first_entry[0]))
        if vcf_dict[key][pos_1][0] == first_entry[0][-1]:
          reverse_flag = True
        if  reverse_flag == True:
          for k in range(len(first_entry)):
            position = int(remove_character(first_entry[k]))
            haplotype_1 = first_entry[k][-1]
            if vcf_dict[key][position][1] == haplotype_1:
              smo.write(str(j)+"\t"+str(1)+"\t"+str(0)+"\t"+str(key)+"\t"+str(position)+"\t"+haplotype_1+"\t"+ vcf_dict[key][position][0]+"\t"+ vcf_dict[key][position][2]+"\n")
              j = j+1
            elif vcf_dict[key][position][0] == haplotype_1:
              smo.write(str(j)+"\t"+str(0)+"\t"+str(1)+"\t"+str(key)+"\t"+str(position)+"\t"+vcf_dict[key][position][1]+"\t"+haplotype_1+"\t"+ vcf_dict[key][position][2]+"\n")
              j = j+1
        else:
          for k in range(len(first_entry)):
            position = int(remove_character(first_entry[k]))
            haplotype_1 = first_entry[k][-1]
            if vcf_dict[key][position][1] == haplotype_1:
              smo.write(str(j)+"\t"+str(0)+"\t"+str(1)+"\t"+str(key)+"\t"+str(position)+"\t"+haplotype_1+"\t"+ vcf_dict[key][position][0]+"\t"+ vcf_dict[key][position][2]+"\n")
              j = j+1
            elif vcf_dict[key][position][0] == haplotype_1:
              smo.write(str(j)+"\t"+str(1)+"\t"+str(0)+"\t"+str(key)+"\t"+str(position)+"\t"+vcf_dict[key][position][1]+"\t"+haplotype_1+"\t"+ vcf_dict[key][position][2]+"\n")
              j = j+1
        smo.write("********\n")
        i= i+1
smo.close()

  
