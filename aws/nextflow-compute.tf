resource "aws_iam_instance_profile" "nf_ecs_instance_role" {
  name = "nextflow-ecs-instance-role"
  role = aws_iam_role.nf_ecs_role.name
}

resource "aws_batch_compute_environment" "nextflow_ondemand" {
  compute_environment_name = "nextflow-ondemand-compute"
  compute_resources {
    instance_role = aws_iam_instance_profile.nf_ecs_instance_role.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "BEST_FIT"
    max_vcpus = 20
    min_vcpus = 0
    image_id = "ami-0a8857ac38c35157f"
    security_group_ids = [
      aws_security_group.nf_security.id,
    ]
    subnets = [
      aws_subnet.nf_subnet.id,
    ]
    type = "EC2"
  }

  service_role = aws_iam_role.nf_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.nf_batch]
}

resource "aws_batch_compute_environment" "nextflow_spot" {
  compute_environment_name = "nextflow-spot-compute"
  compute_resources {
    instance_role = aws_iam_instance_profile.nf_ecs_instance_role.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = aws_iam_role.nf_ecs_role.arn
    bid_percentage = 20
    max_vcpus = 100
    min_vcpus = 0
    image_id = "ami-0a8857ac38c35157f"
    security_group_ids = [
      aws_security_group.nf_security.id,
    ]
    subnets = [
      aws_subnet.nf_subnet.id,
    ]
    type = "EC2"
  }

  service_role = aws_iam_role.nf_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.nf_batch]
}