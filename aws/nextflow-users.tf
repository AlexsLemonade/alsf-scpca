# This file creates a nextflow user and adds it to the nextflow group.
# This user can be used for job submissions, but jobs can also be submitted as individual users.

resource "aws_iam_user" "nf_user" {
  name = "nextflow-batch-user"
  tags = var.default_tags
}

resource "aws_iam_access_key" "nf_key" {
  user = aws_iam_user.nf_user.name
}

resource "aws_iam_user_group_membership" "nf_batch" {
  user = aws_iam_user.nf_user.name
  groups = [
    aws_iam_group.nf_group.name
  ]
}
