## File: makefile
## David E. Taylor 
## Applied Research Laboratory
## Department of Computer Science and Engineering
## Washington University in Saint Louis
## det3@arl.wustl.edu
##
## Synthetic Database Generator 
## See README file for details
 
CC = g++
LIB = -I./ip-address/include
##CFLAGS = -g -pg
CFLAGS = -O2 -std=c++11 $(LIB)
.cc.o:
	${CC} ${CFLAGS} -c $*.cc 

sbintree.o : stdinc.h dlist.h sbintree.h
dbintree.o : stdinc.h dlist.h dbintree.h
random_db.o : stdinc.h FilterList.h random_db.h
custom_db.o : stdinc.h FilterList.h ProtList.h FlagList.h ExtraList.h PortList.h PrefixList.h dlist.h sbintree.h dbintree.h redundant_filter_check.h TupleBST.h custom_db.h ip_prefix.h trie.h filter_graph.h
dlist.o : stdinc.h dlist.h
redundant_filter_check.o : stdinc.h redundant_filter_check.h
## nesting_filter_check.o : stdinc.h nesting_filter_check.h
ProtList.o : stdinc.h ProtList.h
FlagList.o : stdinc.h FlagList.h
ExtraList.o : stdinc.h ExtraList.h
PortList.o : stdinc.h PortList.h
PrefixList.o : stdinc.h PrefixList.h
FilterList.o : stdinc.h FilterList.h
TupleBST.o : stdinc.h dlist.h TupleBST.h
db_generator.o : stdinc.h random_db.h PortList.h FilterList.h custom_db.h
uint128_t.o : uint128_t.h
ip_prefix.o : ip_prefix.h
trie.o : trie.h
flow_network.o : flow_network.h
filter_graph.o : stdinc.h filter_graph.h

db_generator:	db_generator.o random_db.o dlist.o ProtList.o FlagList.o ExtraList.o PortList.o PrefixList.o sbintree.o dbintree.o redundant_filter_check.o FilterList.o TupleBST.o custom_db.o uint128_t.o ip_prefix.o trie.o flow_network.o filter_graph.o
	${CC} ${CFLAGS} db_generator.o random_db.o dlist.o sbintree.o dbintree.o ProtList.o FlagList.o ExtraList.o PortList.o PrefixList.o FilterList.o redundant_filter_check.o TupleBST.o custom_db.o uint128_t.o ip_prefix.o trie.o flow_network.o filter_graph.o -o db_generator

all:	db_generator

clean:
	rm -rf *.o db_generator
