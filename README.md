# diy_examples

Checkout DIY library from : https://github.com/diatomic/diy

This code is tested with the tag 3.5.0.

Update (10/26/2021) : Now the code is adapted for the DIY master branch.

#set diy-include folder

export DIY_INC= (path to diy/include directory)

#Compile:

mpic++ diy_broadcast_merge.cpp -I${DIY_INC} -lpthread

#Run:

mpirun -np 1 ./a.out -d 1 -b 2 -v true
