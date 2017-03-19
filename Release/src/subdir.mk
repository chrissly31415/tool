################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/ToolCalc.cpp \
../src/ToolCalcLJ.cpp \
../src/ToolCalcPP.cpp \
../src/ToolCalcZE.cpp \
../src/ToolIO.cpp \
../src/ToolProp.cpp \
../src/energyFunc.cpp \
../src/toolc.cpp 

OBJS += \
./src/ToolCalc.o \
./src/ToolCalcLJ.o \
./src/ToolCalcPP.o \
./src/ToolCalcZE.o \
./src/ToolIO.o \
./src/ToolProp.o \
./src/energyFunc.o \
./src/toolc.o 

CPP_DEPS += \
./src/ToolCalc.d \
./src/ToolCalcLJ.d \
./src/ToolCalcPP.d \
./src/ToolCalcZE.d \
./src/ToolIO.d \
./src/ToolProp.d \
./src/energyFunc.d \
./src/toolc.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/eigen3 -O2 -funroll-loops -march=native -Wall -Wextra -c -fmessage-length=0 -std=c++11 -fopenmp -pthread -v -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


