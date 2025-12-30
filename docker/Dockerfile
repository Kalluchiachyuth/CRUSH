FROM python:3.10-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /crush

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    pybind11-dev \
    libcurl4-openssl-dev \
    wget \
    curl \
    git \
    bedtools \
    ca-certificates \
    default-jre \
    bc \
    mawk \
    && rm -rf /var/lib/apt/lists/*


RUN pip install --upgrade pip && \
    pip install --no-cache-dir \
        numpy \
        scipy \
        statsmodels \
        tqdm \
        cooler \
        hic-straw

COPY crush-src /crush
RUN chmod +x /crush/crush

ENTRYPOINT ["/crush/crush"]
CMD ["--help"]
