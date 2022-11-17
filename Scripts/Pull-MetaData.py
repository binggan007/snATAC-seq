#Yu Fu
#Weng Lab
#UMass Chan Medical School 
#February 2022

# This script retrives matching bulk ATAC-seq data for single cell ATAC-seq data

import os, json, requests

# This function retrieves biosample information and returns [biosample summary, fileID, simple term name] for each experiment
## takes Experiment accession and T/F bulk assay flag
def getExperimentMetaData(exp,flag):
  url="https://www.encodeproject.org/experiments/"+exp+"/?format=json"
  response = requests.get(url,headers=headers,auth=auth)
  data = response.json()
  biosample_summary = data["simple_biosample_summary"]
  term_name = data["biosample_ontology"]["term_name"]
  ## change if you want different output
  ## return peak for bulk assays 
  if flag:
    for file in data["files"]:
      if file["file_type"] == "bed narrowPeak" and file["assembly"] == genome and "preferred_default" in file and file["preferred_default"] == True:
        return biosample_summary, file["accession"], term_name
  else:
  ## return archr project ID for single cell
    for file in data["files"]:
      if file["file_type"] == "tar" and file["output_type"] == "archr project":
        return biosample_summary, file["accession"], term_name
  return biosample_summary, "---", term_name

# This function store the matching bulk assay experiment accession and file ID for a single cell experiment
## takes the biosample_summary/key of single cell experiments, bulk assay you wanna match, and more parameters.
def matchExperimentBiosample(key,assay,genome,term_name):
  url="https://www.encodeproject.org/search/?type=Experiment" + \
      "&perturbed=false" + \
      "&assay_title=" + assay + \
      "&biosample_ontology.term_name=" + term_name + \
      "&replicates.library.biosample.donor.organism.scientific_name=" + species + \
      "&format=json&limit=all"
  response = requests.get(url, headers=headers,auth=auth)
  data = response.json()

  #search for matching experiment for all biosamples stored in myDict
  for EXP in data["@graph"]:
    biosample_summary,signal,term_name = getExperimentMetaData(EXP["accession"],True)
    suffix = " nuclear fraction" # suffix added to snATAC-seq biosample summary
    phrase = term_name + " " + biosample_summary + suffix
    if phrase==key:
      myDict[key][2] = EXP['accession']
      myDict[key][3] = signal
      return


## Parameters
genome = "GRCh38"
species = "Homo+sapiens"
matchingAssay = "ATAC-seq"
headers = {'accept': 'application/json'}
auth = ('GWZGYYQ3', 'rwv4bpnzsn5ggplx') # personal authentication from ENCODE

## Build query
## only retrieve 20 snATAC-seq datasets as an example
snURL = "https://www.encodeproject.org/search/?type=Experiment" + \
             "&perturbed=false" + \
             "&assay_title=snATAC-seq" + \
             "&replicates.library.biosample.donor.organism.scientific_name=" + species + \
             "&format=json&limit=20" # change to "all" for all experiments
response = requests.get(snURL, headers=headers,auth=auth)
data = response.json()

## print retrieved data
#print(json.dumps(data, indent=4))

print("Pulling "+ str(len(data["@graph"])) + " human snATAC-seq experiments")

## retrieve matching data
myDict = {}
for EXP in data["@graph"]:
  biosample_summary,fileID,term_name = getExperimentMetaData(EXP['accession'],False)
  # structure your dictionary to match the output table you want
  myDict[term_name+" "+biosample_summary] = [term_name,EXP['accession'],"---","---",fileID]
for key in myDict:
  matchExperimentBiosample(key,matchingAssay,genome,myDict[key][0])

## print out the table
print("\t".join(["Biosample","Experiment",matchingAssay+" Experiment",matchingAssay+" Peak","snATAC ArchR Project","Biosample description"]))
for key in myDict:
  print("\t".join(myDict[key])+"\t"+key)
