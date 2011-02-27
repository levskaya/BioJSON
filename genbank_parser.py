import re
import json
#from optparse import OptionParser
from pyparsing import *

#parser = OptionParser()
#parser.add_option("-o", "--output", dest="output",
#                  help="output json file")
#(options, args) = parser.parse_args()

infile=open("danio_rerio_tuba1.gb",'r').read()

locusgrammar = Literal("LOCUS") + \
               Word(alphas+nums+'-_') + \
               Word(nums)+Suppress(Literal('bp')) + \
               Word(alphas+'-') + \
               (Literal("linear")|Literal("circular")).setResultsName("topology") + \
               Suppress(Optional(Word(alphas))) + \
               Word(alphas+nums+'-')    

LPAREN = Suppress("(")
RPAREN = Suppress(")")
SEP = Suppress(Literal(".."))
Integer = Word(nums+"<>").setParseAction(lambda s,l,t: int(t[0].replace("<","").replace(">","")) )
SimpleSlice=Group(Integer + SEP + Integer)

complexSlice = Forward()
complexSlice << (Literal("complement")|Literal("join"))+LPAREN+Group(delimitedList(complexSlice)|delimitedList(SimpleSlice)+RPAREN)
featLocation = Group(SimpleSlice|complexSlice)

Feature = Word(alphas+nums+"_-").setResultsName("type") + \
          featLocation + \
          OneOrMore(Group(Suppress('/') + Word(alphas+nums+"_-") + \
                          Suppress('=') + \
                          (QuotedString('"',multiline=True) | OneOrMore(CharsNotIn("\n")))
                    ))
Features = Suppress(Literal("FEATURES")+Literal("Location/Qualifiers"))+OneOrMore(Feature)

print locusgrammar.searchString(infile)

for feat in Feature.searchString(infile):
    foo=feat.asList()
    fdict = dict(foo[2:])
    fdict["type"]=foo[0]
    fdict["loc"]=foo[1]
    print fdict


Sequence = OneOrMore(Suppress(Word(nums)) + OneOrMore(Word("ACGTacgtNn")))
SequenceRegion = Suppress(Literal("ORIGIN"))+Sequence+Suppress(Literal('//'))

print "".join(SequenceRegion.searchString(infile)[0])





