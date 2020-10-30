# alsf-scpca
Management and analysis tools for ALSF Single-cell Pediatric Cancer Atlas data.

- [Environment Setup](#environment-setup)
  - [Nextflow](#nextflow)
    - [Installing Nextflow](#installing-nextflow)
    - [Nextflow Tower](#nextflow-tower)
  - [Docker](#docker)
  - [AWS](#aws)
    - [Installing AWS Command line tools](#installing-aws-command-line-tools)
    - [Configuring your AWS credentials](#configuring-your-aws-credentials)
    - [AWS infrastructure](#aws-infrastructure)
- [Running Workflows](#running-workflows)
- [Running locally](#running-locally)
  - [Running on AWS Batch](#running-on-aws-batch)

## Environment Setup

### Nextflow

Most of the workflows in this repository are run via [Nextflow](https://www.nextflow.io).

#### Installing Nextflow

You can install Nextflow by following the [instructions on their page](https://www.nextflow.io/docs/latest/getstarted.html#installation).
Be sure to follow the second step of their instructions: move the `nextflow` file to a directory in your `$PATH` for easy access.

Alternatively, you can install `nextflow` via `conda` or `brew` (it may also be available in other package managers):
- `conda install -c bioconda nextflow`
- `brew tap brewsci/bio; brew install nextflow`

#### Nextflow Tower

Optionally, you can track your nextflow workflows in progress with [Nextflow Tower](https://tower.nf/).
To use this, you will need to create an account on https://tower.nf/ (using your github account is likely simplest), then get your token value from https://tower.nf/tokens.
You will then want to create an environment variable with this token value:
`export TOWER_ACCESS_TOKEN=<MYTOWERTOKEN>`
This line can be entered at the command line, but you will likely want to add it to your shell initialization file (`~/.bash_profile`, `~/.zshrc`, etc.)

### Docker

All of the workflows are designed to run via Docker containers for each tool. It is therefore critical that you have Docker installed for any local execution. For downloads and installation, go to https://www.docker.com/products/docker-desktop

### AWS

The tools in this repository are designed primarily to run on AWS via the AWS Batch system, with data stored in S3.
To use these tools as they are currently implemented, you will therefore need an AWS account, with access to the Alex's Lemonade CCDL organization account.

#### Installing AWS Command line tools
Once you have your account information, you will want to install the [AWS command line interface](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html) on your local machine.
You should be able to use either version 1 (>v1.17) or version 2 of the AWS command line tools, but these instructions have been primarily tested with v1.18.
To install, you can follow [amazon's instructions for v2](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html) or install via your favorite package manager (`pip` and `conda` will at this time install v1, `homebrew` will install v2):
- `pip install awscli`
- `conda install -c conda-forge awscli`
- `brew install awscli`

#### Configuring your AWS credentials

When you set up your account, you may have gotten a paired `AWS Access Key ID` and `AWS Secret Access Key`.
If you did not, you can login to the [aws console](https://console.aws.amazon.com/), select "My Security Credentials" from the menu that appears clicking on your username in the upper right.
Scroll down to "Access keys for CLI, SDK, & API access" and click the "Create Access Key" button.
You will see a popup with your Access key ID and a link to show the Secret access key.

In the terminal, run the command `aws configure` and when prompted paste in the `AWS Access Key ID` and `AWS Secret Access Key`.
You may want to set up your `Default region name`: most of the computing resources used here run in `us-east-1`, but setting this to a different region should not affect most commands.
The `Default output format` can be left as `None`

The `aws configure` command will create a directory at `~/.aws` with a `config` file and a `credentials` file that will be required for `nextflow` commands to work smoothly with Batch and S3.

#### AWS infrastructure

The infrastructure on AWS (Batch queues, AMIs, etc.) is defined via [terraform](https://www.terraform.io) using the files in the `aws` directory. More details are described in [aws/setup-log.md](aws/setup-log.md)



## Running Workflows

The workflows for this repository are stored in the `workflows` directory, organized by task into subdirectories.
In general, workflows should be run from the main `workflows` directory, rather than the subdirectories.
This will allow easy access to the shared `nextflow.config` file and keep all intermediate and log files in a single location.
In particular, `nextflow` creates a separate work directory for every task: these will appear by default (for local tasks) in the `workflows/work` directory if the command is invoked from `workflows`.
As the work directories can get large, it is helpful to have a single location to keep track of and purge as needed.

## Running locally

⚠️⚠️⚠️
Running workflows in this repository locally is not something to take lightly!
It is fine for a quick test, and useful for development, but note that most of the workflows here use very large data files, which will have to be downloaded locally to run.
Workflow processing may require large amounts of RAM and time, so the following example commands will rarely be used, and are mostly for illustration.
⚠️⚠️⚠️

The basic command to run a workflow will look something like the following:

```
nextflow run checks/check-md5.nf
```

This will run the workflow locally, using the default parameters as defined in the workflow file.

In most cases, you will want to skip any cached steps that have already run: this can be done by adding the `-resume` flag.

```
nextflow run checks/check-md5.nf -resume
```

⚠️ Again, you probably don't want to run locally unless you have a good reason and know the limitations! ⚠️

### Running on AWS Batch

To run the same workflow on AWS Batch, make sure your credentials are configured [as described above](#configuring-your-aws-credentials), and then run with the `batch` profile, which has been configured in `nextflow.config` for the CCDL infrastructure.

```
nextflow run checks/check-md5.nf -profile batch -resume
```

(Note that the effectiveness `-resume` flag depends on location: locally cached steps will still have to run on AWS, but if they are cached on AWS, they will be skipped.)

If you have set up [Nextflow Tower](#nextflow-tower), you can add a flag to send progress information there:

```
nextflow run checks/check-md5.nf -profile batch -resume -with-tower
```

Finally, if you want to change any of the parameters that are defined in a workflow, you can do so at the command line using flags that start with with a double dash `--`. For example, to use different run ids for the alevin workflow, you might use:

```
nextflow run checks/check-md5.nf -profile batch -resume --run_ids SCPCR000003,SCPCR000004
```