resource "aws_security_group" "nf_security" {
  name = "nextflow-security-group"

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_vpc" "nf_vpc" {
  cidr_block = "10.1.0.0/16"
}

resource "aws_subnet" "nf_subnet" {
  vpc_id = aws_vpc.nf_vpc.id
  cidr_block = "10.1.1.0/24"
}
