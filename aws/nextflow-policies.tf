# S3 Group policies taken from AWS Nextflow batch setup
resource "aws_iam_policy" "nf_readwrite_S3" {
  name   = "nextflow-ccdl-readwrite-s3"
  policy = <<EOF
{
 "Version": "2012-10-17",
 "Statement": [
  {
   "Effect": "Allow",
   "Action": [
    "s3:PutAnalyticsConfiguration",
    "s3:GetObjectVersionTagging",
    "s3:ReplicateObject",
    "s3:GetObjectAcl",
    "s3:GetBucketObjectLockConfiguration",
    "s3:PutLifecycleConfiguration",
    "s3:GetObjectVersionAcl",
    "s3:PutObjectTagging",
    "s3:DeleteObject",
    "s3:DeleteObjectTagging",
    "s3:GetBucketPolicyStatus",
    "s3:GetObjectRetention",
    "s3:GetBucketWebsite",
    "s3:PutReplicationConfiguration",
    "s3:DeleteObjectVersionTagging",
    "s3:PutObjectLegalHold",
    "s3:GetObjectLegalHold",
    "s3:GetBucketNotification",
    "s3:PutBucketCORS",
    "s3:GetReplicationConfiguration",
    "s3:ListMultipartUploadParts",
    "s3:PutObject",
    "s3:GetObject",
    "s3:PutBucketNotification",
    "s3:PutBucketLogging",
    "s3:GetAnalyticsConfiguration",
    "s3:PutBucketObjectLockConfiguration",
    "s3:GetObjectVersionForReplication",
    "s3:GetLifecycleConfiguration",
    "s3:GetInventoryConfiguration",
    "s3:GetBucketTagging",
    "s3:PutAccelerateConfiguration",
    "s3:DeleteObjectVersion",
    "s3:GetBucketLogging",
    "s3:ListBucketVersions",
    "s3:ReplicateTags",
    "s3:RestoreObject",
    "s3:ListBucket",
    "s3:GetAccelerateConfiguration",
    "s3:GetBucketPolicy",
    "s3:PutEncryptionConfiguration",
    "s3:GetEncryptionConfiguration",
    "s3:GetObjectVersionTorrent",
    "s3:AbortMultipartUpload",
    "s3:PutBucketTagging",
    "s3:GetBucketRequestPayment",
    "s3:GetObjectTagging",
    "s3:GetMetricsConfiguration",
    "s3:PutBucketVersioning",
    "s3:GetBucketPublicAccessBlock",
    "s3:ListBucketMultipartUploads",
    "s3:PutMetricsConfiguration",
    "s3:PutObjectVersionTagging",
    "s3:GetBucketVersioning",
    "s3:GetBucketAcl",
    "s3:PutInventoryConfiguration",
    "s3:GetObjectTorrent",
    "s3:PutBucketWebsite",
    "s3:PutBucketRequestPayment",
    "s3:PutObjectRetention",
    "s3:GetBucketCORS",
    "s3:GetBucketLocation",
    "s3:ReplicateDelete",
    "s3:GetObjectVersion"
   ],
   "Resource": [
    "arn:aws:s3:::nextflow-ccdl-data/*",
    "arn:aws:s3:::nextflow-ccdl-results/*",
    "arn:aws:s3:::nextflow-ccdl-data",
    "arn:aws:s3:::nextflow-ccdl-results"
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

resource "aws_iam_policy" "nf_read_S3" {
  name   = "nextflow-ccdl-read-s3"
  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObjectVersionTorrent",
        "s3:GetObjectAcl",
        "s3:GetObject",
        "s3:GetObjectTorrent",
        "s3:GetObjectRetention",
        "s3:GetObjectVersionTagging",
        "s3:GetObjectVersionAcl",
        "s3:GetObjectTagging",
        "s3:GetObjectVersionForReplication",
        "s3:GetObjectLegalHold",
        "s3:GetObjectVersion",
        "s3:ListMultipartUploadParts"
      ],
      "Resource": "arn:aws:s3:::*/*"
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
        "arn:aws:s3:::*",
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