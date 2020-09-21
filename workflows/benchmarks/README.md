## Docker

Running the benchmark analysis notebook should be done via docker, with the following command to mount the local directory as the home direcory and pass in the current user's AWS credentials and launch RStudio.
AWS credentials are expected to be stored in `~/.aws`, having been set up on the user's  system with `aws configure`.
https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html#cli-configure-quickstart-config

```
docker run \
  --mount type=bind,target=/home/rstudio,source=$PWD \
  --mount type=bind,target=/home/rstudio/.aws,source=$HOME/.aws \
  -e PASSWORD=<MYPASS> \
  -p 8787:8787 \
  ghcr.io/alexslemonade/scpca-r
```

The rstudio server can then be logged into via the browser by navigating to http://localhost:8787.
Login with the username `rstudio` and the password chosen above.