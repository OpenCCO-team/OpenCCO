
FROM debian:bullseye

# install additional debian packages
RUN apt-get update && apt-get install -y build-essential cmake libboost-all-dev m4 imagemagick libceres-dev libceres1 git && rm -rf /var/lib/apt/lists/*

# copy the code to $bin
ENV bin /workdir/bin/
RUN mkdir -p $bin
WORKDIR $bin
COPY . .
RUN git clone https://github.com/DGtal-team/DGtal.git; cd DGtal; mkdir build; cd build; cmake .. -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DCMAKE_BUILD_TYPE:string="Release"; make;
RUN rm -rf build; mkdir build &&  cd build; cmake .. -DCMAKE_BUILD_TYPE:string="Release" -DDGtal_DIR="$bin/DGtal/build"  && make generateTree2D generateTree3D mergeObj && cp bin/generateTree2D $bin &&  cp bin/generateTree3D $bin &&  cp tools/mergeObj $bin


# the execution will happen in the folder /workdir/exec
# it will be created by IPOL

# some QoL tweaks
ENV PYTHONDONTWRITEBYTECODE 1
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION python
ENV PATH $bin:$PATH

# $HOME is writable by the user `ipol`, but 
ENV HOME /home/ipol
# chmod 777 so that any user can use the HOME, in case the docker is run with -u 1001:1001
RUN groupadd -g 1000 ipol && useradd -m -u 1000 -g 1000 ipol -d $HOME && chmod -R 777 $HOME
USER ipol