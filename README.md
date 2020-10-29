# alsf-scpca
Management and analysis tools for ALSF Single-cell Pediatric Cancer Atlas data.

- [Environment](#environment)
  - [Installing AWS Command line tools](#installing-aws-command-line-tools)
  - [Configuring your AWS credentials](#configuring-your-aws-credentials)
- [Nextflow](#nextflow)
  - [Installing Nextflow](#installing-nextflow)
  - [Nextflow Tower](#nextflow-tower)

## Environment

The tools in this repository are designed primarily to run on AWS via the AWS Batch system, with data stored in S3.
To use these tools as they are currently implemented, you will therefore need an AWS account, with access to the Alex's Lemonade CCDL organization account.

### Installing AWS Command line tools
Once you have your account information, you will want to install the [AWS command line interface](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html) on your local machine.
You should be able to use either version 1 (>v1.17) or version 2 of the AWS command line tools, but these instructions have been primarily tested with v1.18.
To install, you can follow [amazon's instructions for v2](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html) or install via your favorite package manager (`pip` and `conda` will at this time install v1, `homebrew` will install v2):
- `pip install awscli`
- `conda install -c conda-forge awscli`
- `brew install awscli`

### Configuring your AWS credentials

When you set up your account, you may have gotten a paired `AWS Access Key ID` and `AWS Secret Access Key`.
If you did not, you can login to the [aws console](https://console.aws.amazon.com/), select "My Security Credentials" from the menu that appears clicking on your username in the upper right.
Scroll down to "Access keys for CLI, SDK, & API access" and click the "Create Access Key" button.
You will see a popup with your Access key ID and a link to show the Secret access key.

In the terminal, run the command `aws configure` and when prompted paste in the `AWS Access Key ID` and `AWS Secret Access Key`.
You may want to set up your `Default region name`: most of the computing resources used here run in `us-east-1`, but setting this to a different region should not affect most commands.
The `Default output format` can be left as `None`

The `aws configure` command will create a directory at `~/.aws` with a `config` file and a `credentials` file that will be required for `nextflow` commands to work smoothly with Batch and S3.

## Nextflow

Most of the workflows in this repository are run via [Nextflow](https://www.nextflow.io).

### Installing Nextflow
You can install Nextflow by following the [instructions on their page](https://www.nextflow.io/docs/latest/getstarted.html#installation).
Be sure to follow the second step of their instructions: move the `nextflow` file to a directory in your `$PATH` for easy access.

Alternatively, you can install `nextflow` via `conda` or `brew` (it may also be available in other package managers):
- `conda install -c bioconda nextflow`
- `brew tap brewsci/bio; brew install nextflow`

### Nextflow Tower

Optionally, you can track your nextflow workflows in progress with [Nextflow Tower](https://tower.nf/).
To use this, you will need to create an account on https://tower.nf/ (using your github account is likely simplest), then get your token value from https://tower.nf/tokens.
You will then want to create an environment variable with this token value:
`export TOWER_ACCESS_TOKEN=<MYTOWERTOKEN>`
This line can be entered at the command line, but you will likely want to add it to your shell initialization file (`~/.bash_profile`, `~/.zshrc`, etc.)





