## Docker

Running the benchmark analysis notebook should be done via docker, with the following command to mount the local directory as the home direcory and pass in the current user's AWS credentials and launch RStudio.

```
docker run \
  --mount type=bind,target=/home/rstudio,source=$PWD \
  --mount type=bind,target=/home/rstudio/.aws,source=$HOME/.aws \
  -e PASSWORD=<MYPASS> -p 8787:8787 \
  ccdl/scpca_r
```
