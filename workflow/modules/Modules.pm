package Modules;

# ---------------------------------------------------------------------------- #

# Module name:	Modules
# Created by:	Eliot Stanton
# Created on:	09 September, 2023
# Modified:	13 September, 2023
# Description:	List of modules to be applied before running workflow.

# ---------------------------------------------------------------------------- #

use strict;
use warnings;

use Stats;

use Analysis;
use Check;
use Configuration;
use Comms;
use Edit;
use Forks;
use Stage;
use Wrapper;

use Getopt::Std;
use Parallel::ForkManager;
use Term::ANSIColor qw(:constants);
use YAML::XS;

# ---------------------------------------------------------------------------- #

1;
