This folder contains a Dockerfile for the [vireo-snp](https://vireosnp.readthedocs.io/en/latest/) package.

The image is built from a versioned Python (slim) image.

## Building the image

The image can be built with the following command run in this directory:

```
docker build . -t ghcr.io/alexslemonade/vireo-snp:v0.5.7
```

Alternatively, the latest version should be available via:

```
docker pull ghcr.io/alexslemonade/vireo-snp
```

### Pushing to GitHub Container Registry

If the updated image needs to be pushed to the Github Container Registry, you can follow the steps outlined below (updating version numbers as needed).

Note that you will need to login to the github container registry with a PAT.
See [Github Docs](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry#authenticating-to-the-container-registry) for instructions for logging in.

*Do not push an image unless you are sure it is working!*

The current image was pushed with the following commands.

```
docker tag ghcr.io/alexslemonade/vireo-snp:v0.5.7 ghcr.io/alexslemonade/vireo-snp:latest
docker push ghcr.io/alexslemonade/vireo-snp:v0.5.7
docker push ghcr.io/alexslemonade/vireo-snp:latest
```

