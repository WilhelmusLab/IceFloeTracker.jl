FROM python:3.9

WORKDIR ~/IceFloeTracker.jl/scripts

ENV TERM=xterm

COPY . .

RUN apt-get clean && apt-get update && \
    apt-get install -y wget python3-pip && \
    rm -rf /var/lib/apt/list/*

RUN apt-get install gdal-bin -y

RUN gdalwarp --version

RUN pip3 install skyfield \ 
        requests
RUN echo python --version

RUN bash ./fetchdata.sh -o ./images -s 2022-05-04 -e 2022-05-06 -c wgs84 78.186394 38.250605 70.749318 15.32373