ARG FROM_IMAGE_NAME=nvcr.io/nvidia/openmm:8.1.1

FROM ${FROM_IMAGE_NAME}
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
RUN apt-get update \
    && apt-get install -y git build-essential python3-dev make cmake wget\
    && rm -rf /var/lib/apt/lists/*
    
# RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh -O /tmp/miniconda.sh
# RUN /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
#     rm /tmp/miniconda.sh && \
#     echo "export PATH=/opt/conda/bin:$PATH" > /etc/profile.d/conda.sh
# ENV PATH /opt/conda/bin:$PATH

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/anaconda.sh
RUN /bin/bash ~/anaconda.sh -b -p "/opt/conda"
RUN . /opt/conda/etc/profile.d/conda.sh && \
    . /opt/conda/bin/activate
ENV PATH="/opt/conda/bin:$PATH:"
RUN conda update --all && \
    conda install -y -c conda-forge pdbfixer openmmforcefields biopython
ADD . .
