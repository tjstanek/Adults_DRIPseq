## Alternate approach:  Start with empty dictionaries, use if/else statements to add new/append counts to corresponding features
# Import test list of shuffled reads
import sys
fh = open(sys.argv[1])
x = fh.readline()
y = x.rstrip()
z = y.split()
feature = z[0]
fh.close()

fh = open(sys.argv[1])

# Set counters
featcounts = 0

# Create empty dictionaries for Autosome counts and X chrom counts
Alldict = {}

for line in fh:
    line = line.rstrip()
    lis = line.split()
    if lis[0] == feature: # Match line feature name to feature variable
        featcounts += int(lis[2]) # Append values to Autosome counter
        
    else: # When feature name changes, append counter values to dictionary keys before changing features and resetting counters
        
        if feature in Alldict: # find key matching current feature
            Alldict[feature].append(featcounts) # append Autosome counter to list of matching key
        else: # if key is not yet present in dictionary
            Alldict[feature] = [featcounts] # Create new key-value pair with Autosome counter as a list
        featcounts = int(lis[2])
    feature = lis[0] # Update feature with new feature name

if feature in Alldict: # find key matching current feature
    Alldict[feature].append(featcounts) # append Autosome counter to list of matching key
else: # if key is not yet present in dictionary
    Alldict[feature] = featcounts # Create new key-value pair with Autosome counter as a list

fh.close()


fh = open(sys.argv[2])
outfile = open(sys.argv[3],"w")
Allcounts = {}
for line in fh:
    line = line.rstrip()
    lis = line.split()
    Allcounts[lis[0]] = int(lis[2])

for k,v, in Alldict.items():
    gtcount = 0
    if k in Allcounts:
        vlis = list(v)
        for x in vlis:
            if x > Allcounts[k]:
                gtcount += 1
        outfile.write(k + "\tAll\t" + str(gtcount/10000) + "\n")
        gtcount = 0

fh.close()
outfile.close()
