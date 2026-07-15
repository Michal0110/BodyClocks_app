FROM rocker/r-ver:4.5.2

ENV DEBIAN_FRONTEND=noninteractive
ENV RENV_CONFIG_PAK_ENABLED=FALSE

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    pkg-config \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libicu-dev \
    libglpk-dev \
    zlib1g-dev \
    libuv1-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy dependency metadata first, allowing Docker to cache package installation.
COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json
COPY .Rprofile .Rprofile

RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

RUN R -e "renv::restore(prompt = FALSE)"

# Copy the application only after dependencies have been installed.
COPY . .

EXPOSE 7860

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=7860)"]