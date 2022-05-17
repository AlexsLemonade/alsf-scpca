resource "aws_launch_template" "nf_lt_standard" {
  name = "nextflow-launchtemplate-standard"
  tags = var.default_tags
  # ccdl-nextflow-base-v2.0 image
  image_id = "ami-01c08d4de548477df"
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size = 128
      encrypted = true
      delete_on_termination = true
    }
  }
  update_default_version = true
}

resource "aws_launch_template" "nf_lt_bigdisk" {
  name = "nextflow-launchtemplate-bigdisk"
  tags = var.default_tags
  # ccdl-nextflow-base-v2.0 image
  image_id = "ami-01c08d4de548477df"
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size = 1000
      encrypted = true
      delete_on_termination = true
    }
  }
  update_default_version = true
}

resource "aws_launch_template" "nf_lt_auto_scaled_ebs" {
  name = "nextflow-launchtemplate-auto-scaled-ebs"
  tags = var.default_tags
  # ccdl-nextflow-base-v2.0 image
  image_id = "ami-01c08d4de548477df"
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size           = 40  # GiB (image size + 10 GiB)
      encrypted             = true
      delete_on_termination = true
    }
  }
  update_default_version = true
  user_data = filebase64("./user_data/auto_scaled_ebs.yaml")
}
