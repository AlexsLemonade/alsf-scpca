This folder contains a Dockerfile for the scpca-aws image used for scripted interactions with S3.

The image is built from a versioned miniconda3 image, with boto3 and the aws-cli tools installed.

The image can be built with the following command run in this directory:

```
docker build . -t ghcr.io/alexslemonade/scpca-aws
```

Alternatively, the latest version should be available via:

```
docker pull ghcr.io/alexslemonade/scpca-aws
```

