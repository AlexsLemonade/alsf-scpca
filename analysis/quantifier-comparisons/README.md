
## Docker setup

Running the quantification comparison notebooks are most easily done via Docker, using the `ghcr.io/alexslemonade/scpca-r` image, which contains and RStudio server and all necessary dependencies, including the AWS command line tools and R packages.  AWS credentials are expected to be stored in `~/.aws` on the user's system, having been set up with `aws configure` on the host. https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html#cli-configure-quickstart-config

Use the following command to mount a local directory as the home directory within the image (it is recommended that you use the repository root: `alsf-scpca`), pass in the current user's AWS credentials, and launch RStudio server (substitute a password for <MYPASS>).

```
cd alsf-scpca
docker run \
  --mount type=bind,target=/home/rstudio,source=$PWD \
  --mount type=bind,target=/home/rstudio/.aws,source=$HOME/.aws \
  -e PASSWORD=<MYPASS> \
  -p 8787:8787 \
  ghcr.io/alexslemonade/scpca-r
```

The RStudio server can then be logged into via the browser by navigating to http://localhost:8787.
Login with the username `rstudio` and the password chosen above.

Periodically the latest version of the Docker image may need to be refreshed with:
```
docker pull ghcr.io/alexslemonade/scpca-r
```