FROM continuumio/miniconda3:latest

RUN mkdir /opt/taranys/
ADD utils /opt/taranys/utils
ADD test /opt/taranys/test
ADD *.py /opt/taranys/
ADD environment.yml /opt/taranys/
ADD logging_config.ini /opt/taranys/
ADD README.md /opt/taranys/
ADD LICENSE /opt/taranys/

SHELL ["/bin/bash", "-c"]
RUN cd /opt/taranys
RUN /opt/conda/bin/conda env create -f /opt/taranys/environment.yml && /opt/conda/bin/conda clean -a
RUN /opt/conda/bin/conda env export --name taranys > taranys.yml
RUN echo "conda activate taranys" > ~/.bashrc
ENV PATH /opt/conda/envs/taranys:/opt/conda/envs/taranys/utils:$PATH
