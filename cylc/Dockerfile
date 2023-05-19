FROM julia:1.9

RUN apt-get clean && apt-get update && \
    apt-get install -y wget python3-pip && \
    rm -rf /var/lib/apt/list/*

RUN apt-get install gdal-bin -y

RUN gdalwarp --version

# RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.0-linux-x86_64.tar.gz
# RUN tar zxvf julia-1.9.0-linux-x86_64.tar.gz
# RUN cp -r julia-1.9.0 /opt/
# RUN ln -s /opt/julia-1.9.0/bin/julia /usr/local/bin/julia
RUN julia -e 'using Pkg'

RUN pip3 install skyfield \ 
        requests
RUN python --version

# RUN julia activate ./scripts
# COPY ./scripts