This folder contains a Dockerfile for the [alevin-fry](https://github.com/COMBINE-lab/alevin-fry) image.

The image is built from a versioned miniconda3 image, with Alevin-fry obtained from bioconda and additional utilities needed to run Nextflow. 
 
You can build the image with the following command:

```
docker build . -t ghcr.io/alexslemonade/scpca-alevin-fry
```

Alternatively, the latest version should be available via: 
```
docker pull ghcr.io/alexslemonade/scpca-alevin-fry
```
