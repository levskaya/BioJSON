"""
GenBank FlatFile Parser

Parses GenBank FlatFile (GFF) into a JSON representation that's quickly usable in other scripts / etc.

Copyright (c) 2011 Anselm Levskaya (http://anselmlevskaya.com)
Licensed under the MIT (http://www.opensource.org/licenses/mit-license.php) license.  
"""
import re
import json
from optparse import OptionParser
from pyparsing import *

#===============================================================================
# GenBank Grammar
#===============================================================================

# Generic Entry
#===============================================================================

# All entries in a genbank file headed by an all-caps title with no space between start-of-line and title
CapWord = Word("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
# after titled line, all subsequent lines have to have at least one space in front of them
# this is how we split up the genbank record
SpacedLine =  White(min=1) + CharsNotIn("\n") + LineEnd()
#HeaderLine = CapWord + CharsNotIn("\n") + LineEnd()
GenericEntry =  Dict(Group(CapWord + Combine(CharsNotIn("\n") + LineEnd() +\
                             ZeroOrMore( SpacedLine ))))


# GenBank LOCUS Entry Parser
#===============================================================================
# LOCUS string, the first line of any genbank-sh file
# unlike many shite parsers, this should work for NCBI, ApE, and NTI style gb files
LocusEntry =   Literal("LOCUS") + \
               Word(alphas+nums+'-_').setResultsName("name") + \
               Word(nums)+Suppress(Literal('bp')) + \
               Word(alphas+'-').setResultsName("moleculetype") + \
               (Literal("linear")|Literal("circular")).setResultsName("topology") + \
               Suppress(Optional(Word(alphas))) + \
               Word(alphas+nums+'-').setResultsName("date")

# GenBank Feature Table Parser
#===============================================================================

#==== Genbank Location String Parser
#
# a string of slices w. functional modifiers that go at most two levels deep
# slice is N1..N2  w. N1<N2
# i.e.
# 23..88  --> seq[23:89] in python syntax (genbank uses inclusive slicing)
# 234..555 
# complement(234..434) --> rc(seq[234:435])
# join(23..343,454..666,777..999) --> seq[23:344]+seq[454:667]+seq[777:1000]
# complement(join(23..343,454..666,777..999))
# join(complement(34..123),complement(333..565))
#
# additionally the slices can have ambiguous locs like <454..999 or 232..>331
# also note the dumb 34.38 fuzzy slice notation
# i.e. <45..900  says the real feature starts "somewhere" before pos 45
#       45.48 says feature somewhere between bases 45->48
# lot of weird annotations best avoided because the connote ~useless knowledge for synthetic design
#
# if you don' know where something is, don't use it or guess and move on

LPAREN = Suppress("(")
RPAREN = Suppress(")")
SEP = Suppress(Literal(".."))

#recognize numbers w. < & > uncertainty specs, then strip the <> chars to make it fixed
gbIndex = Word(nums+"<>").setParseAction(lambda s,l,t: int(t[0].replace("<","").replace(">","")) )
SimpleSlice=Group(gbIndex + SEP + gbIndex)

complexSlice = Forward()
complexSlice << (Literal("complement") | Literal("join")) + LPAREN + ( delimitedList(complexSlice) | delimitedList(SimpleSlice) ) + RPAREN 
featLocation = Group( SimpleSlice | complexSlice)

def parseGBLoc(s,l,t):
    """retwingles parsed genbank location strings,
    assumes no joins of RC and FWD sequences """
    strand = 1
    locationlist=[]
    
    #see if there are any complement operators
    for entry in t[0]:
        if entry=="complement": strand=-1

    for entry in t[0]:
        if type(entry)!=type("string"):
            locationlist.append([entry[0],entry[1],strand])
            
    #return locationlist
    return [['location', locationlist ]]

featLocation.setParseAction(parseGBLoc)

#==== Genbank Feature Key-Value Pairs

# NCBI feature format:   /key="value"
NCBIFeaturekeyval = Group(Suppress('/') + Word(alphas+nums+"_-") + Suppress('=') + QuotedString('"',multiline=True) )                    
# ApE doesn't put its vals in quotes!   /key=value  (I'm assuming it doesn't do multilines this way?)
ApEFeaturekeyval = Group(Suppress('/') + Word(alphas+nums+"_-") + Suppress('=') + OneOrMore(CharsNotIn("\n")) )

Feature = Group( Word(alphas+nums+"_-").setResultsName("type").setParseAction(lambda s,l,t: [ ["type", t[0]] ] ) +\
                 featLocation.setResultsName("location") +\
                 OneOrMore( NCBIFeaturekeyval | ApEFeaturekeyval ) )

FeaturesEntry = Literal("FEATURES") + Literal("Location/Qualifiers") + Group(OneOrMore(Feature)).setResultsName("features")


# GenBank Sequence Parser
#===============================================================================
# sequence is just a column-spaced big table of dna nucleotides
# should recognize full IUPAC alphabet
Sequence = OneOrMore(Suppress(Word(nums)) + OneOrMore(Word("ACGTacgtNn")))

# Group(  ) hides the setResultsName names def'd inside, such that one needs to first access this group and then access the dict of contents inside
SequenceEntry = Group(Literal("ORIGIN") + Sequence.setResultsName("sequence").setParseAction(lambda s,l,t: "".join(t) )).setResultsName("sequence")


# Final GenBank Parser
#===============================================================================

#GB files with multiple records split by "//" sequence at beginning of line
GBEnd = Literal("//")

GB = LocusEntry + OneOrMore(FeaturesEntry | SequenceEntry | GenericEntry) + GBEnd

#===============================================================================
# End Genbank Parser
#===============================================================================

if __name__ == "__main__":
    #parse command line string 
    parser = OptionParser()
    #parser.add_option("-o", "--output", dest="output",
    #                  help="output json file")
    (options, args) = parser.parse_args()
    
    infile = open(args[0],'r').read()

    parsed = GB.parseString(infile)
    
    jseq = { "name" : parsed["name"],
             "type" : parsed["moleculetype"],
             "date" : parsed["date"],
             "topology" : parsed["topology"],
             "sequence" : parsed["sequence"]["sequence"],
             "features" : map(dict,parsed['features'].asList())
             }


    print json.dumps(jseq)

    outfile=open(args[1],'w')
    json.dump(jseq,outfile)


# Cruft that needs moved to unit tests
"""
for k in ("name","topology","date","moleculetype"):
    print k + ": " + parsed[k]

for f in parsed['features']:
    print "-"*20
    df=dict(f.asList())
    for k in df.keys():
        print k, " : ", df[k]

 print featLocation.parseString("complement(join(34..99,555..999))")
 print featLocation.parseString("join(complement(34..99),complement(555..999))")
 print featLocation.parseString("34..55")
 print featLocation.parseString("complement(34..55)")
 print featLocation.parseString("join(34..99,555..999,1234..7777)")
 print parsed
"""
