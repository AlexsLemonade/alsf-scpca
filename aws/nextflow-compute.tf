# This file creates the compute environments used by the default and priority queues
# The default environment is a 100 vCPU spot cluster
# Priority environment is a 20 vCPU on demand cluster

resource "aws_iam_instance_profile" "nf_ecs_instance_role" {
  name = "nextflow-ecs-instance-role"
  tags = var.default_tags
  role = aws_iam_role.nf_ecs_role.name
}

# Create an spot instance environment with up to 256 vcpus
# the AMI used is described in setup-log.md
resource "aws_batch_compute_environment" "nf_spot" {
  compute_environment_name = "nextflow-spot-compute"
  tags = var.default_tags
  compute_resources {
    instance_role = aws_iam_instance_profile.nf_ecs_instance_role.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = aws_iam_role.nf_spotfleet_role.arn
    bid_percentage = 100
    max_vcpus = 256
    min_vcpus = 0
    # standard launch template
    launch_template {
      launch_template_id = aws_launch_template.nf_lt_standard.id
      version = aws_launch_template.nf_lt_standard.latest_version
    }
    # ec2_key_pair = aws_key_pair.nf_keypair.key_name
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
        parent = "nextflow-spot-compute"
      }
    )
  }

  service_role = aws_iam_role.nf_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.nf_batch_role]
}

# Create an spot instance0 environment with up to 128 vcpus with large disks

resource "aws_batch_compute_environment" "nf_spot_bigdisk" {
  compute_environment_name = "nextflow-spot-compute-bigdisk"
  tags = var.default_tags
  compute_resources {
    instance_role = aws_iam_instance_profile.nf_ecs_instance_role.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = aws_iam_role.nf_spotfleet_role.arn
    bid_percentage = 100
    max_vcpus = 128
    min_vcpus = 0
    # large disk launch template
    launch_template {
      launch_template_id = aws_launch_template.nf_lt_bigdisk.id
      version = aws_launch_template.nf_lt_bigdisk.latest_version
    }
    # ec2_key_pair = aws_key_pair.nf_keypair.key_name
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
        parent = "nextflow-spot-compute"
      }
    )
  }

  service_role = aws_iam_role.nf_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.nf_batch_role]
}

# Create a spot instance environment with up to 2048 vcpus and auto-scaled EBS.
resource "aws_batch_compute_environment" "nf_spot_auto_scaled_ebs" {
  compute_environment_name = "nextflow-spot-compute-auto-scaled-ebs"
  tags = var.default_tags
  compute_resources {
    instance_role = aws_iam_instance_profile.nf_ecs_instance_role.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = aws_iam_role.nf_spotfleet_role.arn
    bid_percentage = 100
    max_vcpus = 2048
    min_vcpus = 0

    launch_template {
      launch_template_id = aws_launch_template.nf_lt_auto_scaled_ebs.id
      version = aws_launch_template.nf_lt_auto_scaled_ebs.latest_version
    }
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
        parent = "nextflow-spot-compute"
      }
    )
  }

  service_role = aws_iam_role.nf_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.nf_batch_role]
}

# Create an ondemand environment with up to 512 vcpus
resource "aws_batch_compute_environment" "nf_ondemand" {
  compute_environment_name = "nextflow-ondemand-compute"
  tags = var.default_tags
  compute_resources {
    instance_role = aws_iam_instance_profile.nf_ecs_instance_role.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "BEST_FIT"
    max_vcpus = 512
    min_vcpus = 0
    # standard launch template
    launch_template {
      launch_template_id = aws_launch_template.nf_lt_auto_scaled_ebs.id
      version = aws_launch_template.nf_lt_auto_scaled_ebs.latest_version
    }
    # ec2_key_pair = aws_key_pair.nf_keypair.key_name
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
        parent = "nextflow-ondemand-compute"
      }
    )
  }
  service_role = aws_iam_role.nf_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.nf_batch_role]
}
