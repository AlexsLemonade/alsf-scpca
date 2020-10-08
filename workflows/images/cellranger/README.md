This folder contains a Dockerfile for the cellranger analysis.

The cellranger archive component must be downloaded separately to comply with licensing, and should be placed in this folder.
The current version of this file is `cellranger-4.0.0.tar.gz`


```
docker build . -t scpca-cellranger:4.0.0
```

This is then tagged and pushed to AWS ECR (must login to ECR first):
```
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 589864003899.dkr.ecr.us-east-1.amazonaws.com
docker tag scpca-cellranger:4.0.0 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
docker tag scpca-cellranger:4.0.0 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:4.0.0
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:4.0.0
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest


```
docker pull 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
```

Note that this is a private repository; access only to Alexslemonade users at this time.