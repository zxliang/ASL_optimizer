CXX		:= g++
CFLAGS		:= -std=c++11

SRC_DIR		:= src
SRCS		:= $(shell find $(SRC_DIR) -name '*.cpp')

OBJ_DIR		:= build
OBJS		:= $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

TARGET		:= bin/ASL_optimizer
INC		:= -I ./include

#all:
#	@echo $(SRCS)
#	@echo $(OBJS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(INC) -c $< $(CFLAGS) -o $@

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(CFLAGS)

.PHONY : clean
clean:
	rm -f $(TARGET) $(OBJS)



