make config MFEM_CXX="nvcc" \
CXXFLAGS="-Xcompiler -fopenmp -O3 --restrict --expt-extended-lambda -x=cu -arch=sm_70 -std=c++11 -m64" \
MFEM_EXT_LIBS="-L/usr/local/cuda/lib64-lrt -lcuda -lcudart -lcudadevrt -lnvToolsExt -lgomp -L/usr/WS1/vargas45/COUPLED_2019/fix_proxy_dev/RAJA/build/lib -lRAJA" \
MFEM_TPLFLAGS=-I/usr/WS1/vargas45/COUPLED_2019/fix_proxy_dev/RAJA/build/include \
OPTIM_FLAGS=-O3 MFEM_USE_MM=YES MFEM_USE_RAJA=YES
