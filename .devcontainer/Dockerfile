FROM ubuntu:latest

RUN apt-get update && apt-get install -y vim gcc g++ gdb make cmake gfortran libblas-dev liblapack-dev pkgconf git wget python3 python3-pip

# Install fpm
# WORKDIR /tmp
# RUN git clone --depth=1 https://github.com/fortran-lang/fpm
# WORKDIR /tmp/fpm
# RUN wget https://github.com/fortran-lang/fpm/releases/download/current/fpm.F90
# RUN mkdir -p build/bootstrap
# RUN gfortran -J build/bootstrap -o build/bootstrap/fpm fpm.F90
# RUN ./build/bootstrap/fpm install --prefix=/usr/local
# WORKDIR /root

# Install stdlib
# RUN fpm install --prefix /usr/local stdlib

# Install fypp
WORKDIR /tmp
RUN git clone --depth=1 https://github.com/aradi/fypp.git
WORKDIR /tmp/fypp
RUN cp ./src/fypp.py /usr/local/bin/fypp
WORKDIR /root

# Install stdlib
WORKDIR /tmp
RUN git clone --depth=1 https://github.com/fortran-lang/stdlib
WORKDIR /tmp/stdlib
RUN mkdir build
WORKDIR /tmp/stdlib/build
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=off ..
RUN make install -j 8
WORKDIR /root
