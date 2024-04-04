# Specific policies used by roles and groups
# Which S3 buckets are available for reading

# S3 Group policies taken from AWS Nextflow batch setup

# This policy allows read and write access to specific buckets for nextflow processing
resource "aws_iam_policy" "nf_readwrite_S3" {
  name   = "nextflow-ccdl-readwrite-s3"
  tags   = var.default_tags
  policy = <<EOF
{
 "Version": "2012-10-17",
 "Statement": [
  {
   "Effect": "Allow",
   "Action": [
    "s3:ListBucket",
    "s3:*Object"
   ],
   "Resource": [
    "arn:aws:s3:::nextflow-ccdl-data",
    "arn:aws:s3:::nextflow-ccdl-data/*",
    "arn:aws:s3:::nextflow-ccdl-results",
    "arn:aws:s3:::nextflow-ccdl-results/*",
    "arn:aws:s3:::openscpca-temp-simdata",
    "arn:aws:s3:::openscpca-temp-simdata/*"
   ]
  },
  {
   "Effect": "Allow",
   "Action": [
    "s3:GetAccountPublicAccessBlock",
    "s3:ListAllMyBuckets",
    "s3:ListAccessPoints",
    "s3:HeadBucket"
   ],
   "Resource": "*"
  }
 ]
}
EOF
}

# This policy gives read access to all S3 buckets
resource "aws_iam_policy" "nf_read_S3" {
  name   = "nextflow-ccdl-read-s3"
  tags   = var.default_tags
  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject*",
        "s3:ListBucket"
      ],
      "Resource": [
        "arn:aws:s3:::ccdl-scpca-data",
        "arn:aws:s3:::ccdl-scpca-data/*",
        "arn:aws:s3:::nextflow-ccdl-data",
        "arn:aws:s3:::nextflow-ccdl-data/*",
        "arn:aws:s3:::nextflow-ccdl-results",
        "arn:aws:s3:::nextflow-ccdl-results/*"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetLifecycleConfiguration",
        "s3:GetBucketTagging",
        "s3:GetInventoryConfiguration",
        "s3:ListBucketVersions",
        "s3:GetBucketLogging",
        "s3:ListBucket",
        "s3:GetAccelerateConfiguration",
        "s3:GetBucketPolicy",
        "s3:GetEncryptionConfiguration",
        "s3:GetBucketObjectLockConfiguration",
        "s3:GetBucketRequestPayment",
        "s3:GetAccessPointPolicyStatus",
        "s3:GetMetricsConfiguration",
        "s3:GetBucketPublicAccessBlock",
        "s3:GetBucketPolicyStatus",
        "s3:ListBucketMultipartUploads",
        "s3:GetBucketWebsite",
        "s3:GetBucketVersioning",
        "s3:GetBucketAcl",
        "s3:GetBucketNotification",
        "s3:GetReplicationConfiguration",
        "s3:DescribeJob",
        "s3:GetBucketCORS",
        "s3:GetAnalyticsConfiguration",
        "s3:GetBucketLocation",
        "s3:GetAccessPointPolicy"
      ],
      "Resource": [
        "arn:aws:s3:::ccdl-scpca-data",
        "arn:aws:s3:::nextflow-ccdl-data",
        "arn:aws:s3:::nextflow-ccdl-results",
        "arn:aws:s3:*:*:accesspoint/*",
        "arn:aws:s3:*:*:job/*"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetAccessPoint",
        "s3:GetAccountPublicAccessBlock",
        "s3:ListAllMyBuckets",
        "s3:ListAccessPoints",
        "s3:ListJobs",
        "s3:HeadBucket"
      ],
      "Resource": "*"
    }
  ]
}
EOF
}

resource "aws_iam_policy" "nf_manage_ebs" {
  name        = "nextflow-ccdl-manage-ebs"
  description = "A policy that allows to manage (attach/create/delete) EBS volumes."
  tags        = var.default_tags
  policy      = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "ec2:AttachVolume",
                "ec2:DescribeVolumeStatus",
                "ec2:DescribeVolumes",
                "ec2:DescribeTags",
                "ec2:ModifyInstanceAttribute",
                "ec2:DescribeVolumeAttribute",
                "ec2:CreateVolume",
                "ec2:DeleteVolume",
                "ec2:CreateTags"
            ],
            "Resource": "*"
        }
    ]
}
EOF
}
