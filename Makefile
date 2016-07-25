FILES = $(wildcard *.cpp)
CC = g++
COMMANDS = -o
OUTPUT_NAME = testoutput

OUTPUT_NAME : $(FILES)
	$(CC) $(FILES) $(COMMANDS) $(OUTPUT_NAME)
