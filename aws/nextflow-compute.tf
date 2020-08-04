resource "aws_iam_instance_profile" "nf_ecs_instance_role" {
  name = "nextflow-ecs-instance-role"
  role = aws_iam_role.nf_ecs_role.name
}

# Create an ondemand environment with up to 20 vcpus
# the AMI used is described in setup-log.md
resource "aws_batch_compute_environment" "nf_ondemand" {
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
    tags = merge(
      var.default_tags,
      {
        parent = "nextflow-batch-ondemand"
      }
    )
  }

  service_role = aws_iam_role.nf_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.nf_batch]

}

# Create an spot instance environment with up to 100 vcpus
# the AMI used is described in setup-log.md
resource "aws_batch_compute_environment" "nf_spot" {
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
    type = "SPOT"
    tags = merge(
      var.default_tags,
      {
        parent = "nextflow-batch-ondemand"
      }
    )
  }

  service_role = aws_iam_role.nf_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.nf_batch]
}
