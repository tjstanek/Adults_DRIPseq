import sys
print("Command line arguments are", sys.argv[1], "and", sys.argv[2])
#Set initial feature, to be modified as code iterates through file
fh = open(sys.argv[1])
x = fh.readline()
y = x.rstrip()
z = y.split()
feature = z[0]
fh.close()

# Re-open file to loop through lines and assign counts to features
fh = open(sys.argv[1])
# Set counter
featcounts = 0

# Create empty dictionaries for Autosome counts and X chrom counts
Alldict = {}

#Set initial feature, to be modified as code iterates through file
for line in fh:
    line = line.rstrip()
    lis = line.split()
    if lis[0] == feature: # Match line feature name to feature variable
        featcounts += int(lis[2]) # Append values to Autosome counter
        
    else: # When feature name changes, append counter values to dictionary keys before changing features and resetting counters

        if feature in Alldict: # find key matching current feature
            #print(feature) # sanity check that key/feature match is working
            Alldict[feature].append(featcounts) # append Autosome counter to list of matching key
        else: # if key is not yet present in dictionary
            #print(feature) # sanity check that new features are being added
            Alldict[feature] = featcounts # Create new key-value pair with Autosome counter as a list
        featcounts = int(lis[2])
    feature = lis[0] # Update feature with new feature name

if feature in Alldict: # find key matching current feature
    #print(feature) # sanity check that key/feature match is working
    Alldict[feature].append(featcounts) # append Autosome counter to list of matching key
else: # if key is not yet present in dictionary
    #print(feature) # sanity check that new features are being added
    Alldict[feature] = featcounts # Create new key-value pair with Autosome counter as a list

fh.close()

# Combine all output counts  of desired *.cf and shuffle.*.cf files using features template above

final=open(sys.argv[2],"w")
with open("DRIPseq_all_features_overlap.tab") as fh2:
      for feat in fh2: # open DRIPseq features
        feat=feat.rstrip()
        fis = feat.split()
        if fis[0] in Alldict.keys():
            pass
        else:
            Alldict[fis[0]] = str(0)
for k,v, in Alldict.items():
    final.write(k + '\tAll\t' + str(v) + '\n')

final.close()
