This folder contains a Dockerfile for the scpca-r image used for R analysis in this repository.

The image is built from a versioned Rocker image, with additional R packages and other utilities added as needed.

It can be built with the following command run in this directory:

```
docker build . -t ghcr.io/alexslemonade/scpca-r
```

Alternatively, the lastest version should be available via:

```
docker pull ghcr.io/alexslemonade/scpca-r
```

To use and launch the RStudio server, you can use the following example command, which uses the current directory for the internal home directory.
Note that this also mounts the current user's `.aws` credentials folder, which can be useful for easy S3 access, as most data files in this repository are stored on S3 in private buckets.

```
docker run \
  --mount type=bind,target=/home/rstudio,source=$PWD \
  --mount type=bind,target=/home/rstudio/.aws,source=$HOME/.aws \
  -e PASSWORD=<YOURPASS> \
  -p 8787:8787 \
  ghcr.io/alexslemonade/scpca-r
```