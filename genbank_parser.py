import re
import json
from optparse import OptionParser
from pyparsing import *

parser = OptionParser()
parser.add_option("-o", "--output", dest="output",
                  help="output json file")
(options, args) = parser.parse_args()

infile=open(args[0],'r').read()

#GenBank Grammar

# All entries in a genbank file headed by an all-caps title with no space between start-of-line and title
CapWord = Word("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
# after titled line, all subsequent lines have to have at least one space in front of them
# this is how we split up the genbank record
SpacedLine =  White(min=1) + CharsNotIn("\n") + LineEnd()
#HeaderLine = CapWord + CharsNotIn("\n") + LineEnd()
GenericEntry =  CapWord + Combine(CharsNotIn("\n") + LineEnd() +\
                             ZeroOrMore( SpacedLine ))
GBEnd = Literal("//")
#GB = OneOrMore(GenericEntry) + GBEnd

# LOCUS string, the first line of any genbank-sh file
# unlike many shite parsers, this should work for NCBI, ApE, and NTI style gb files
LocusEntry =   Literal("LOCUS") + \
               Word(alphas+nums+'-_') + \
               Word(nums)+Suppress(Literal('bp')) + \
               Word(alphas+'-') + \
               (Literal("linear")|Literal("circular")).setResultsName("topology") + \
               Suppress(Optional(Word(alphas))) + \
               Word(alphas+nums+'-')    

# Feature Location
# a string of slices w. functional modifiers that go at most two levels deep
# slice is N1..N2  w. N1<N2
# i.e.
# 234..555
# complement(234..434)
# join(23..343,454..666,777..999)
# complement(join(23..343,454..666,777..999))
# join(complement(34..123),complement(333..565))
# ...additionally the slices can have ambiguous locs like <454..999 or 232..>331 also dumb 34.38 fuzzy slice
# lot of weird annotations best avoided because the connote ~useless knowledge for synthetic design
LPAREN = Suppress("(")
RPAREN = Suppress(")")
SEP = Suppress(Literal(".."))
gbIndex = Word(nums+"<>").setParseAction(lambda s,l,t: int(t[0].replace("<","").replace(">","")) )
SimpleSlice=Group(gbIndex + SEP + gbIndex)
complexSlice = Forward()
complexSlice << (Literal("complement") | Literal("join")) + LPAREN + Group( delimitedList(complexSlice) | delimitedList(SimpleSlice) ) + RPAREN 
featLocation = Group( SimpleSlice | complexSlice)
#feature definition
Feature = Word(alphas+nums+"_-") + featLocation + \
          OneOrMore(Group(Suppress('/') + Word(alphas+nums+"_-") + \
                          Suppress('=') + \
                          # ApE doesn't put its vals in quotes...
                          (QuotedString('"',multiline=True) | OneOrMore(CharsNotIn("\n")))
                    ))

FeaturesEntry = Literal("FEATURES") + Literal("Location/Qualifiers") + Group(OneOrMore(Feature)).setResultsName("features")

Sequence = OneOrMore(Suppress(Word(nums)) + OneOrMore(Word("ACGTacgtNn")))
SequenceEntry = Group(Literal("ORIGIN") + Sequence.setResultsName("sequence").setParseAction(lambda s,l,t: "".join(t) )).setResultsName("se")

GB = LocusEntry + OneOrMore(FeaturesEntry | SequenceEntry | GenericEntry) + GBEnd

print "total parse test"
parsed = GB.parseString(infile)
print parsed
print
print parsed['features']
print parsed['se']['sequence']

#print locusgrammar.searchString(infile)

# for feat in Feature.searchString(infile):
#     foo=feat.asList()
#     fdict = dict(foo[2:])
#     fdict["type"]=foo[0]
#     fdict["loc"]=foo[1]
#     print fdict
#     print

#print "all entries"
#print dict(Splitgram.searchString(infile).asList())
#print CapWord.searchString(infile)
#print Spaceleaded.searchString(infile)


