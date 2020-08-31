#!/usr/bin/env python
#
# Function to convert list of genotypes into list of possible site patterns
#  for D(FOIL) analysis
# Heterozygosity can exist in any sample except P0
# No missing data allowed!

# topology: (((P1,P2),(P3,P4)),P0)

P0gt = "0/0"
P1gt = "1/0"
P2gt = "1/0"
P3gt = "1/0"
P4gt = "0/0"

def GTs2sitepattern (P0, P1, P2, P3, P4):
    GTs = [P1,P2,P3,P4,P0]
    ### get 'A' state from P0
    A = P0.split("/")[0]
    ### remove redundant alleles
    for i in range(5):
        GT = GTs[i].split("/")
        if GT[0] == GT[1]:
            if GT[0] == A:
                GTs[i] = ["A"]
            else:
                GTs[i] = ["B"]
        else:
            GTs[i] = ["A","B"]
    ### LOOP THROUGH GENOTYPES TO PRODUCE A LIST OF SITE PATTERNS
    patterns = []
    for a in range(len(GTs[0])):
        for b in range(len(GTs[1])):
            for c in range(len(GTs[2])):
                for d in range(len(GTs[3])):
                    patterns.append("".join([GTs[0][a],
                                  GTs[1][b],
                                  GTs[2][c],
                                  GTs[3][d]])+"A")
    return(patterns)

print(GTs2sitepattern(P0gt, P1gt, P2gt, P3gt, P4gt))
