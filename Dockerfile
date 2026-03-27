FROM condaforge/mambaforge:24.9.0-0

WORKDIR /opt/all_in_one_pipeline
COPY . /opt/all_in_one_pipeline

RUN mamba env create -y -f env/all_in_one_requirements.yml && \
    conda clean -a -y

ENV PATH="/opt/conda/envs/all_in_one_pipeline/bin:/opt/conda/bin:${PATH}"
ENV CHECKM_DATA_PATH=/export/database/checkm

CMD ["bash"]

