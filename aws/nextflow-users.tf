resource "aws_iam_user" "nextflow" {
  name = "nextflow-batch-user"
}

resource "aws_iam_access_key" "nextflow" {
  user = aws_iam_user.nextflow.name
}

resource "aws_iam_user_group_membership" "nextflow-batch" {
  user = aws_iam_user.nextflow.name
  groups = [
    "${aws_iam_group.nextflow.name}"
  ]
}