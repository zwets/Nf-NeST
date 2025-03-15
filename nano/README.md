# Nanopore branch of nfNeST

This is just a quick hack to make the workflow run with Nanopore reads.

It has not been validated, and probably could be much improved by using
minimap2 or some other Nanopore-optimised mapper.

In particular, a likely much better approach is that adapt EPI2ME's
[wf-tb-amr workflow](https://github.com/epi2me-labs/wf-tb-amr).

All I've done in this repository is copy the original workflow definition,
and modify the steps up to the alignment to take single fastq files.

I have changed no parameters in any of the steps apart from the bbduk
trimming, where I reduced the minimal read Q-score (from 35 to 20) as
otherwise hardly any reads would be included.

Running is quite simple, see the `run-workflow.sh` script in this directory.

