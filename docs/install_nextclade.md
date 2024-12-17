# Nextclade CLI

## Installing Nextclade CLI

[Nextclade](https://clades.nextstrain.org/) is an open-source tool that performs several tasks including sequence alignment, translation, mutation calling, phylogenetic placement, and more. It has a web app version and a CLI tool that is run locally on the computer. 

VARGRAM uses [Nextclade CLI](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/index.html) to create the mutation profile from input sequences, capturing Nextclade's [analysis output file](https://docs.nextstrain.org/projects/nextclade/en/stable/user/output-files/04-results-tsv.html). 

You can use VARGRAM without installing Nextclade CLI. Simply upload your sequences to [Nextclade Web]((https://clades.nextstrain.org/)) and download the analysis output file. Then, feed this analysis file as an input to VARGRAM. However, it is recommended to install Nextclade CLI, which would simplify your workflow.

To install Nextclade CLI, run the necessary command from [the Nextclade tutorial](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/installation/standalone.html) based on your operating system. The commands download the latest Nextclade CLI version from GitHub, renames the Nextclade file to `nextclade` and makes it executable.

Finally, Nextclade has to be added to the path.

## Adding an executable to the path

The path is a list of directories of executable files. To add the executable file `nextclade` to the path, navigate to your terminal emulator (Terminal on Mac and Command Prompt on Windows). On macOS or Linux, you can view the path by typing `echo $PATH` on the command line. On Windows, type `echo %PATH%` instead. Then, follow the instructions below based on your operating system. Once your Nextclade executable is added to the path, Nextclade CLI commands can be run anywhere, which is what VARGRAM assumes.

=== "Linux &  macOS"

    1. Determine whether you are using `zsh` or `bash`. Type `echo $SHELL` to know. 
    2. Run `echo 'export PATH="$PATH:/path/to/your/nextclade/directory/"' >> ~/.zshrc` if you are using `zsh` or `echo 'export PATH="$PATH:/path/to/your/nextclade/directory/"' >> ~/.bashrc` if you are using `bash`.
    3. Run `source ~/.zshrc` for `zsh` or `source ~/.bashrc` for `bash`.
    4. To verify, run `echo $PATH` and you should see the directory of your Nextclade executable. 

=== "Windows"

    1. Open the Settings app.
    2. Navigate to `System > About > Advanced system settings`.
    3. Click on `Environment Variables`.
    4. In the "System variables" section, select `Path` and click on `Edit`.
    5. Click on `New` and then type or paste the path of the directory containing the Nextclade executable.
    6. Click `OK` to close the pop-up windows.
    7. To verify, open the Command Prompt and run `echo %PATH%` and you should see the directory of the Nextclade executable.

!!! warning "Security prompt"

    Your computer may prevent Nextclade from running as a security measure. On macOS, for example, you may get a warning that `nextclade` is from an unidentified developer. Simply grant the necessary permission to `nextclade` on your computer's settings app or on the security prompt that pops up.
