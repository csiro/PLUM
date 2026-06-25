# Code assist using Claude in Amazon Web Services
created by Williams, Gareth

The aim is for everybody to do something hands-on. Setup will be required.

An admin (probably Williams, Gareth) will need to add you to at least one cloud tenancy as noted in: SC Cloud and Coding Assist Evaluation Project

You will almost certainly need somewhere to run commands (or you can use the vendor web-gui and consoles). Your local/windows/mac/linux setup may be fine, but for commonality and some sand-boxing, DiSP is preferred. We can use a terminal in vscode, jupyter, Rstudio or the new test cli-workbench.

## DiSP

DiSP - go to https://rms.it.csiro.au and set up a new project (or use an existing one) and be able to make kubeflow workspaces.

* DiSP vscode - prefer custom image:  
latest preview:  
https://docker-registry.it.csiro.au/disp-runai-test/vscode-uv-ubuntu24-wbtool  
https://docker-registry.it.csiro.au/dkp-images/disp-vscode:beta-v8.1.7

* DiSP cli-workbench - choose jupyter with custom image: 
https://docker-registry.it.csiro.au/scce/util/cli-dev-workbench:v1.11 
* DiSP Rstudio - works but needs more environment tweaks

## Extra CLI tools

Once you have a workspace, from a terminal clone the project at: csiro-internal/rootless-tools: simplify install of a set of tools that do not require root (https://github.com/csiro-internal/rootless-tools)

Note. it would be nice if git setup were more trivial.

If you can't git clone, we can put it somewhere you can scp/rsync.
Today you can (may need user@):

```yaml
scp -rp petrichor.hpc.csiro.au:/tmp/rootless-tools ./
```

You can run the rootless-tools/bin/wb-tool script to install extra tools (as a non-privileged user, into .local/bin).

You probably want 'aws' and 'claude' (and to follow the wb-tool claude tips). If you install kubectl (or already have it) you can potentially run new k8s things too. 

Install:

./rootless-tools/bin/wb-tool install claude
./rootless-tools/bin/wb-tool install aws

Set env vars as recommended by Claude code install script:


```yaml
export AWS_REGION=ap-southeast-2
export CLAUDE_CODE_USE_BEDROCK=1
```


One easy way of working with AWS and the available cloud tenancy is to go to: 

[AWS Tenancy Address](https://csiro.awsapps.com/start/#/console?account_id=811936640717) 

then select 

```yaml
button.name("PowerUserAccess")
``` 

to get the AWS web gui/console with a SSO authenticated session.

Connect to that AWS session in the cli environment with):
```yaml
bash> aws login --remote

and as prompted open the url provided, find the session and copy/paste the token (Shift-Insert is your friend for pasting).
```

Open Claude (first time runs setup):
```yaml
bash> claude
```