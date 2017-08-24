#!/bin/bash

# Installation script for HGTector

set -eu

if ! hash R 2>/dev/null
then
    echo "R is not detected. You may either install R manually and run this script again, or simply omit it."
    exit 1
fi

export PERL_MM_USE_DEFAULT=1

if ! perl -MApp::cpanminus -e 1 2>/dev/null
then
    echo "Installing cpanminus, for easy installation of Perl modules..."
    perl -MCPAN -e 'install App::cpanminus'
    echo done
fi

if ! perl -MStatistics::R -e 1 2>/dev/null
then
    echo "Installing Perl module Statistics::R, for crosstalking between Perl and R..."
    cpanm Statistics::R
    echo done
fi

if ! perl -MSpreadsheet::WriteExcel -e 1 2>/dev/null
then
    echo "Installing Perl module Spreadsheet::WriteExcel, for generating Excel spreadsheets..."
    cpanm Spreadsheet::WriteExcel
    echo done
fi

if [[ $(R -q -e "require('pastecs')" 2>&1) == *"there is no package called"* ]]
then
    echo "Installing R package pastecs, for calculating local extrema..."
    R -e "install.packages('pastecs', repos='https://cran.rstudio.com/')"
    echo done
fi

if [[ $(R -q -e "require('diptest')" 2>&1) == *"there is no package called"* ]]
then
    echo "Installing R package diptest, for performing Hartigan's dip test for unimodality..."
    R -e "install.packages('diptest', repos='https://cran.rstudio.com/')"
    echo done
fi
