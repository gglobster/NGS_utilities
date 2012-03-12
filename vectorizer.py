# script to represent sequences as vectors (based on gene annots)



from sets.nheA_set import ann_30K as genomes


symbolDB = []
vectorDB = []

for genome in genomes:

    g_vector = []

    # load genome file to extract features (to proteins in mfas file)

    # first go
    if len(symbolDB) < 1:
        print "A"
        # add all features as new symbols

    # normal
    else:
        print "B"
        # blast each feature

        # for each feature:

            # if there is a hit above threshold:
                # add feature to symbolDB
                # add hit symbol to vector

            # else:
                # create new symbol entry in symbolDB from feature
                # add new symbol to vector






