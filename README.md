# tristan-mp-PU_densdecomp
Modified version of `TRISTAN-MP` code.  Main features:

   - Esirkepov current deposit for 1st, 2nd, and 3rd order shape functions, in
   2D and 3D
   - Ghost cell modifications, to accommodate higher-order particle shapes

A few summary documents are included in `/doc`:

   - Esirkepov_change_v2s.pdf: document to summarize differences between the modified
   code, and the original (see `/doc/tristan-mp-pu-master.zip`)
   - Esirkepov_current_deposit_v2.pdf: document to summarize several tests of the
   modified code


# tristan-mp-PU
TRISTAN-MP parallel electromagnetic 3D particle-in-cell code.
Pre-release version for testing. 

Developed by: Anatoly Spitkovsky, Luis Gargate, Jaehong Park, Lorenzo Sironi. 
Based on original TRISTAN code by Oscar Buneman. 

See http://tristan-mp.wikispaces.com for more extensive documentation. 
Warning: some of the info on the wiki is obsolete. 

Organization and suggested workflow: 
The code is modular and is organized in such a way that in most cases 
the user would not need to edit any of the main source files in the main 
code directory. 

All of the user-specific configuration should be possible to confine
to user_* routines in "user_yourname" subdirectory. There are example 
user configuration files in default "user" directory, showing 
a counterstreaming Weibel instability setup, two-stream instability 
and a collisionless shock simulation. 
There are also sample input files.

When using code from github, clone it to your local machine, 
switch to "master" branch. You can create your branch off of master.
We recommend making your own user_yourname directory with your specific files.
There will be some Princeton specific files in user_princeton directory. 

On local machine:
git clone https://github.com/PrincetonUniversity/tristan-mp-pu.git

copy some example files to start with

setup your branch:
> git checkout -b yourname

> mkdir user_yourname

> cp user/user_weibel.F90 user_yourname/user_mysetup.F90

edit user_yourname/user_mysetup.F90

To compile:

From the main directory:

> make USER_FILE=user_yourname/user_mysetup

(no need to specify the extension F90)
This command sets the environment variable USER_FILE which 
is the name of your configuration file. You normally do not need 
to edit the Makefile. 

to clean: 
> make clean 

The executable tristan-mp2d will appear in exec/ directory

3D version is enabled when 3D=yes flag is added to the make line, e.g.,
make USER_FILE=... 3D=yes 
 
You need to have parallel HDF5 library installed with intel or GNU compilers, 
which will create h5pfc alias for the compiler. Some instructions for 
installation are on wiki page. For Macs brew seems to work fine with gfortran:

$ brew install gcc
$ brew install openmpi --enable-mpi-thread-multiple
$ brew install hdf5 --with-fortran --with-mpi

To run:
Make a run directory somewhere outside the git-controlled source directory. 
Copy the executable tristan-mp2d there. 
Copy example submit and input files from directory 
(see wiki page for example submit
files for clusters; you don't need these on your desktop/laptop).
 
Input file has to be named "input" in the run directory, or the executable takes -i option. 
E.g.: 

./tristan-mp2d -i input.weibel

(for MPI, it can be, e.g.: srun -n 16 ./tristan-mp2d -i input.weibel)
Note that you need to edit the input file to set the number of domain sub-partitions
 sizex * sizey * sizez be equal to the total number of cores to be used. sizez = 1 in 2D. 

It is possible to select the output directory with -o flag, e.g.:
./tristan-mp2d -i input.weibel -o output.weibel 

Edit submit and input files for your case and according to the queue policy of your cluster. 
In 3D, the domain is split in y and z directions, 
with "sizey" chunks in y direction and "total # of cpus/sizey" chunks in z direction. 

>qsub submit 
or other appropriate submission command. 

When running, the code will create subdirectories output and restart.
The output is single HDF5 files (single per time step). 
Currently we provide routines to interactively view output using python.
https://github.com/pcrumley/Iseult

It requires anaconda to run. We had good experience with anaconda 4.0.0 
on Mac, but not later. The older version is available on anaconda's website. 

To launch the vis tool, run (path to Iseult)Iseult/iseult.py .
This will open interactive windows. Right click on plots to get more options. 

load.py is a script that loads HDF5 files into a python dictionary, 
that can be accessed as d[i]['bz'], where i is the file number. 

There are also older IDL routines, which are available on request. 

Suggested workflow with git and on github:

Normally, you will work in your branch. Branches are selectable
by 
> git checkout yourname

(presuming that "yourname" is the name of your branch)

On the first invocation of your branch you will need to create it:
> git checkout -b yourname master

This will create branch "yourname" based off of master, and will switch 
to it. Then create your directory
> mkdir user_yourname

and populate it with a user file. 

You will frequently record your changes in commits by: 
(from the main directory)

> git add user_yourname

Most of the time all your changes are in this directory, and this 
should be sufficient. You can check git status to see if anything
else important was changed. There is .gitignore that is set up to 
avoid all .o and executable files.  

Usually you won't need more than these additional commands: 
> git add code

> git add user_pu

Of course, you can also do 
> git add -A

if really necessary
 
Commit your changes like this:
> git commit -m "Descriptive message of the change"

Periodically, pull in changes from the master with:
> git pull origin master

or 

> git rebase origin master

(this replays your changes after changes from master are pulled). 

Now, let's say you have made contributions to code or user_pu directory
in your branch that would be useful to others. This is the flow to submit
a pull request for modifying the master:

> git checkout master

> git pull origin master 

> git checkout -b yourname-pull-request 

> git diff master yourname > my-change-set.diff

(I'm assuming "yourname" is the name of your branch 
with the relevant changes)

> cd code/ 

> git apply ../my-change-set.diff

(repeat this for other directories affected, such as user_pu, or user as needed)

> cd <repo-root>

> rm my-change-set.diff

> git add .

> git commit -m "message" 

> git push origin yourname-pull-request

Why could we not just "git merge yourname" into this branch? 
Because this would bring in your user_yourname directory
into master. Removing it later may have adverse consequences 
when you do next "git pull master" on your branch. 

Now go to github.com/PrincetonUniversity/tristan-mp-pu and click on Pull Request, 
selecting the right branch (yourname-pull-request). Others will review your code
and incorporate it into the master.

When all is done and your pull request is approved you can delete this branch:
> git branch -d yourname-pull-request 

More info on pull requests is here:
https://yangsu.github.io/pull-request-tutorial/

Good luck! 


