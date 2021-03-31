This folder contains a Dockerfile for the alevin-fry analysis.


## Building the image

 The image is built from a versioned miniconda3 image, with Alevin-fry obtained from bioconda and additional utilities needed to run Nextflow. 
 
 You can build the image running the following command from this `images/alevinfry` working directory:

```
docker build . -t scpca-alevin-fry:0.1.0
```

At this point, the image should be ready for use on the local machine.

## Using the AWS ECR

To use the AWS ECR docker repository, you will need to login  with:
```
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 589864003899.dkr.ecr.us-east-1.amazonaws.com
```

Note that for the code above to work, you must have set up `aws` command line tools on your machine, and have run [`aws configure`](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html) to set up your credentials for AWS access.
This is a private repository; access is only available to Alex's Lemonade (CCDL) users at this time.

### Pulling the current image

To pull the current version of the container (skipping the build steps), you can use the following command (once you are logged in):

```
docker pull 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-alevin-fry:latest
```

### Pushing to AWS ECR

If the updated image needs to be pushed to AWS ECR, you can follow the outline steps below (updating version numbers as needed).
*Do not push an image unless you are sure it is working!*

The current image was pushed with the following commands.

```
docker tag scpca-alevin-fry:0.1.0 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-alevin-fry:latest
docker tag scpca-alevin-fry:0.1.0 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-alevin-fry:0.1.0
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-alevin-fry:0.1.0
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-alevin-fry:latest
```