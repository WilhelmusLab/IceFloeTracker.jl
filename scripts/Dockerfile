FROM python:3.9-bullseye

ENV TERM=xterm

RUN apt-get clean && apt-get update && \
    apt-get install -y wget python3-pip && \
    rm -rf /var/lib/apt/list/*

RUN apt-get install gdal-bin -y

RUN pip3 install skyfield \ 
        requests

COPY fetchdata.sh /usr/local/bin
COPY pass_time_snakemake.py /usr/local/bin

ENTRYPOINT bash