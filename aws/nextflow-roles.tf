
### Batch Role
resource "aws_iam_role" "nf_batch_role" {
  name = "nextflow-batch-service-role"
  assume_role_policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
    {
        "Action": "sts:AssumeRole",
        "Effect": "Allow",
        "Principal": {
        "Service": "batch.amazonaws.com"
        }
    }
    ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "nf_batch" {
  role = aws_iam_role.nf_batch_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
}


### ECS Role
resource "aws_iam_role" "nf_ecs_role" {
  name = "nextflow-ecs-instance-role"
  assume_role_policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
    {
        "Action": "sts:AssumeRole",
        "Effect": "Allow",
        "Principal": {
        "Service": "ec2.amazonaws.com"
        }
    }
    ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "ecs_ec2_container" {
  role = aws_iam_role.nf_ecs_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_iam_role_policy_attachment" "ecs_rw_s3" {
  role = aws_iam_role.nf_ecs_role.name
  policy_arn = aws_iam_policy.nf_readwrite_S3.arn
}

resource "aws_iam_role_policy_attachment" "ecs_read_s3" {
  role = aws_iam_role.nf_ecs_role.name
  policy_arn = aws_iam_policy.nf_read_S3.arn
}


### Spotfleet Role
resource "aws_iam_role" "nf_spotfleet_role" {
  name = "nextflow-spotfleet-role"
  assume_role_policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
    {
        "Action": "sts:AssumeRole",
        "Effect": "Allow",
        "Principal": {
        "Service": "spotfleet.amazonaws.com"
        }
    }
    ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "nf_spotfleet_tagging" {
  role = aws_iam_role.nf_spotfleet_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole"
}

resource "aws_iam_role_policy_attachment" "nf_spotfleet_autoscale" {
  role = aws_iam_role.nf_spotfleet_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetAutoscaleRole"
}

