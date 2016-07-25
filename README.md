# HaloPlex Cpipe Configuration

## Overview

HaloPlex data has special characteristics that require custom tools and configuration
steps in comparison to a conventional exome or targeted sequencing analysis. This 
repository contains source code for components and pipeline changes needed to customise 
[Cpipe](http://cpipeline.org) for use with HaloPlex data. 


## Cpipe Customisations

The customisations include:

 * Add custom trimming of FASTQ files to remove adapter contamination
 * Remove PCR deduplication step from Cpipe
 * Disable some downstream checks for metrics that are not applicable to Cpipe 
   (eg: duplication rate)


## Applying the Customisations 

To run the analysis Cpipe, first install Cpipe by following directions at 

http://cpipeline.org

Then add the HaloPlex settings in cpipe/designs/DSD3/DSD3.settings.txt to an
analysis profile for the Cpipe installation. Samples configured using the 
DSD3 analysis profile will be analysed using the above customised steps.



