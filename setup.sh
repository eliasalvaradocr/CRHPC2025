#!/bin/bash

module purge

module load cuda/11.8
module load gcc/11.1.0
module load cmake/3.26.0
module load kokkos/cuda

if [ ! -e build/ ]; then
  mkdir build/
fi

CC=gcc CXX=g++ cmake -S . -B build/ -Wno-dev
cmake --build build/ -j8

is_important_file ()
{
  case $1 in
        CMakeFiles)
          ;;
        cmake_install.cmake)
          ;;
        Makefile)
          ;;
        run.slurm)
          return 42
          ;;
        *)
          return 43
          ;;
  esac
}

DIRS=`ls build/basics/`

for DIR in $DIRS
do
  COMPLETE_PATH=build/basics/$DIR
  is_important_file $DIR
  RESULT=$?
  if [ $RESULT -eq 43 -a -d $COMPLETE_PATH ]; then
    FILES=`ls $COMPLETE_PATH`
    for FILE in $FILES
    do
      is_important_file $FILE
      RESULT=$?
      if [ $RESULT -eq 43 ]; then
        FULLPATH=$COMPLETE_PATH/$FILE
        if [ ! -d $FULLPATH ]; then
          echo Copying run.slurm to $PWD/$COMPLETE_PATH ...
          cp basics/run.slurm $COMPLETE_PATH
          echo $FULLPATH
          sed -i "3 s/<program>/${FILE}/" $COMPLETE_PATH/run.slurm
          echo $PWD/$FULLPATH >> $COMPLETE_PATH/run.slurm
        else
          EXAMPLE_FILES=`ls $COMPLETE_PATH/$FILE`
          echo Copying run.slurm to $COMPLETE_PATH/$FILE ...
          cp basics/run.slurm ${COMPLETE_PATH}/${FILE}
          for EXAMPLE_FILE in $EXAMPLE_FILES
          do
            is_important_file $EXAMPLE_FILE
            RESULT=$?
            if [ $RESULT -eq 43 ]; then
              case $EXAMPLE_FILE in
                heat_simulation*)
                  sed -i "3 s/<program>/${EXAMPLE_FILE}/" $COMPLETE_PATH/$FILE/run.slurm
                  echo $PWD/$COMPLETE_PATH/$FILE/$EXAMPLE_FILE 10 10 1000 10 >> $COMPLETE_PATH/$FILE/run.slurm
                  ;;
                kmeans*)
                  sed -i "3 s/<program>/${EXAMPLE_FILE}/" $COMPLETE_PATH/$FILE/run.slurm
                  echo $PWD/$COMPLETE_PATH/$FILE/$EXAMPLE_FILE $PWD/basics/k_means/dataset.csv 4 10 >> $COMPLETE_PATH/$FILE/run.slurm
                  ;;
                Spmv*)
                  sed -i "3 s/<program>/${EXAMPLE_FILE}/" $COMPLETE_PATH/$FILE/run.slurm
                  echo $PWD/$COMPLETE_PATH/$FILE/$EXAMPLE_FILE $PWD/basics/Spvm/data.txt >> $COMPLETE_PATH/$FILE/run.slurm
                  ;;
                *)
                  echo $EXAMPLE_FILE
                  sed -i "s/<program>/${EXAMPLE_FILE}/" $COMPLETE_PATH/$FILE/run.slurm
                  ;;
              esac
            fi
          done
        fi
      fi
    done
  fi
done
