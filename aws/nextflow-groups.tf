resource "aws_iam_group" "nextflow" {
  name = "nextflow-batch"
}

# Batch access
resource "aws_iam_group_policy_attachment" "batch_access" {
  group = aws_iam_group.nextflow.name
  policy_arn = "arn:aws:iam::aws:policy/AWSBatchFullAccess"
}

# EC2 access (may not be needed?)
# resource "aws_iam_group_policy_attachment" "ec2_access" {
#   group      = aws_iam_group.nextflow.name
#   policy_arn = "arn:aws:iam::aws:policy/AmazonEC2FullAccess"
# }

resource "aws_iam_group_policy_attachment" "read_s3" {
  group = aws_iam_group.nextflow.name
  policy_arn = aws_iam_policy.nextflow_read_S3.arn
}

resource "aws_iam_group_policy_attachment" "rw_s3" {
  group = aws_iam_group.nextflow.name
  policy_arn = aws_iam_policy.nextflow_readwrite_S3.arn
}