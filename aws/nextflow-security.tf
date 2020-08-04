resource "aws_security_group" "nf_security" {
  name = "nextflow-security-group"
  vpc_id = aws_vpc.nf_vpc.id
  tags = var.default_tags

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  ingress {
    description = "SSH from anywhere."
    from_port   = 22
    to_port     = 22
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_key_pair" "nf_keypair" {
  key_name = "nextflow-key"
  public_key = "PUT_YOUR_PUBLIC_KEY_HERE"
}

resource "aws_vpc" "nf_vpc" {
  cidr_block = "10.1.0.0/16"
  tags = var.default_tags
}

resource "aws_subnet" "nf_subnet" {
  vpc_id = aws_vpc.nf_vpc.id
  cidr_block = "10.1.1.0/24"
  tags = var.default_tags
}
