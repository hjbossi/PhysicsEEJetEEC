
default: all

all: prepare library/Messenger.o library/BasicUtilities.o library/TauHelperFunctions3.o library/CATree.o library/Dictionary.o library/DrawRandom.o

prepare:
	mkdir -p library/

library/BasicUtilities.o: source/BasicUtilities.cpp include/BasicUtilities.h
	g++ source/BasicUtilities.cpp -Iinclude -c -fpic -o library/BasicUtilities.o -I${RootMacrosBase}/ -std=c++11

library/TauHelperFunctions3.o: source/TauHelperFunctions3.cpp include/TauHelperFunctions3.h
	g++ source/TauHelperFunctions3.cpp -Iinclude -c -o library/TauHelperFunctions3.o -I${RootMacrosBase}/ -std=c++11

library/CATree.o: source/CATree.cpp include/CATree.h
	g++ source/CATree.cpp -Iinclude -c -fpic -o library/CATree.o -I${RootMacrosBase}/ -std=c++11

library/DrawRandom.o: source/DrawRandom.cpp include/DrawRandom.h
	g++ source/DrawRandom.cpp -Iinclude -c -fpic -o library/DrawRandom.o -I${RootMacrosBase}/ -std=c++11

library/Dictionary.o: include/Dictionary.h include/DictionaryObject.h
	rootcint -f source/Dictionary.cxx -c -Iinclude DictionaryObject.h Dictionary.h
	g++ `root-config --cflags` -Iinclude -c source/Dictionary.cxx -o library/Dictionary.o

library/Messenger.o: source/Messenger.cpp include/Messenger.h
	g++ source/Messenger.cpp -Iinclude -c -o library/Messenger.o `root-config --cflags`

