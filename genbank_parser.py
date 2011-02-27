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

#gbtest ="LOCUS       EU482404                 598 bp    mRNA    linear   INV 04-MAR-2008"
#apetest="LOCUS                  New_DNA        2963 bp ds-DNA   linear       13-JAN-2011"
#ntitest="LOCUS       F-GFAP-gCaMP3-W       11512 bp    DNA     circular     14-JUN-2010"
#print locusgrammar.parseString(apetest)
#print locusgrammar.parseString(gbtest)
#print locusgrammar.parseString(ntitest)

SimpleSlice=Group(Word(nums+"<>") + \
                  Suppress(Literal("..")) + \
                  Word(nums+"<>")
                  )

LPAREN = Suppress("(")
RPAREN = Suppress(")")
complexSlice = Forward()
complexSlice << (Literal("complement")|Literal("join"))+LPAREN+Group(delimitedList(complexSlice)|delimitedList(SimpleSlice)+RPAREN)
featLocation = Group(SimpleSlice|complexSlice)

#print featLocation.parseString("345..6565")
#print featLocation.parseString("complement(345..6565)")
#print featLocation.parseString("join(345..6565,7888..9999,12340..13000)")
#print featLocation.parseString("join(complement(12..56),complement(99..1039))")
#print featLocation.parseString("complement(join(12..56,99..1039))")

Feature = Word(alphas+nums+"_-").setResultsName("type") + \
          featLocation + \
          OneOrMore(Group(Suppress('/') + Word(alphas+nums+"_-") + \
                          Suppress('=') + \
                          (QuotedString('"',multiline=True) | OneOrMore(CharsNotIn("\n")))
                    ))
Features = Suppress(Literal("FEATURES")+Literal("Location/Qualifiers"))+OneOrMore(Feature)


print Features.parseString(apefeatures)




