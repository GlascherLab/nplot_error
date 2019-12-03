# nplot_error

NPLOT_ERROR plots data with error bars in various formats and offers the
choice betweeen within-subject and between-subject errorbars

USAGE: [hdl,cfg] = nplot_error('option','value',...)

nplot_error take option value pairs as inputs and returns the axis handle
for further manual plot commands and the config struct (cfg), which can
be also edited manually and used as input for the next nplot_error call.

INPUTS
OPTION           (POSSIBLE) VALUES
'config'         struct with config options (from previous nplot_error)

'type'           'bar','dot','line','fact'

'data'           nSubjects*nConditions data matrix
'mean'           vector of (manually) computed means
'error'          vector of (manually) computed errors

'errorclass'     'bw' (between-subject) or 'ws' (within-subject)
'errortype'      'sem','sd','ci' (90%), 'none'
'color'          cell array of color specs (1*nCondition) [bar,dot,line] OR
                 vector of dimension 1*nCondition with values between
                 1:8 referring to the spcialized color palette defined in
                 cfg.colorder (1=blue, 2=orange, 3=red, 4=green, 5=violet,
                 6=yellow, 7=magenta, 8=brown)
                 Colors can be specified using standard MATLAB characters
                 (rgbcmykw) or as RGB vectors with values between 0 and 1
'errorcolor'     string or numerical color spec for errorbars
                 special string 's' for "same as bar/dot"

'indivdata'      plot individual data points as little dots ('on','off')
'subjectnum'     plot subject numbers next to dots for individual data
                 ('on','off')

'group'          cell array of vectors of subject numbers (rows in data
                 matrix), specify xpos, xtick, and xlabel for one group
                 only, it will be automatically replicated for the remaining
                 groups
'gnames'         cell array of names for groups

'factor'         vector of number of level per factor, left-most factor
                 rotates slowest [type=fact,prod(factor)==size(data,2)]
'factornames'    cell array of factor names [type=fact]
'levelnames'     cell array of cell array of levelnames (for each factor)
                 [type=fact]

'xpos'           positions on x-axis for plotting the data, use this for
                 visual grouping
'xtick'          ticks on x-axis for adding labels for data
'xlabel'         cell array of labels for xtick [type=bar/dot/line]
                 (must correspond to length(xtick)
'xaxislabel'     label for x-axis (type=bar/dot/line]
'yaxislabel'     label for y-axis
'title'          plot title

'legend'         flag to turn on legend ('on','off')
'legnames'       cell array for different colors [type=bar/dot/line]
'leglocation'    location for legend (e.g. 'nw','ne','se','sw')

'figcmd'         command for figure/axes/subplots

nplot_error offers a new plot type: 'fact'orial, which parcellates a data
matrix according to a factor level specification ('factor'). The
left-most factor rotates slowest. The right-most factor rotates fastest and
is plotted on the x-axis. Other factors are shown as different color or
symbols and linestyles, which are automatically assigned. This factor
specification must be valid for the data matrix specified in 'data'.
Factorial plots automatically use a within-subject error. Ticklabels,
axis labels, and legend labels are automatically extracted from the
specifications in 'factornames' and 'levelnames'.

If a within-subject error is also needed for one of the other plot types,
then it needs to be explicitly selected ('errorclass','ws'). In addition,
the 'factor' option has to be specified by usually supplying a single
number ['factor',size(data,2)].

'xpos' can be used for visual grouping along the x-axis, e.g.
'xpos',[1:3 5:7] leaves a small gap between the first 3 and the second 3
conditions. If the option 'group' is specified, then entire plot is
"replicated" for each group along the x-axis. The chosen errorclass
(between- or within-subjects) is used for each of the groups.
_________________________________________________________________________
Within-subject errors are computed according to Morey (2008) [3]. Data
are first normalized w.r.t. to the subject's mean and then standard SEMs
are computed on the normalized data. These SEM/SD/CI are corrected for
the number of factor levels according to Morey (2008). While these
errorbars are not identical to Loftus & Masson (2004) [1], as detailed in
Franz & Loftus (2012) [3], they do provide a good approximation, which allows
for different errorbars for each condition (the Loftus & Masson errorbars
use the pooled error term and hence errorbar are identical for all
levels). The method of paired differences (Franz & Loftus, 2012) [4]
leaves some ambiguities, which errorbar should be used for a condition as
all conditions occur in several paired differences. For a full discussion
see
[1] Loftus & Masson (1994). Using confidence intervals in within-subject
    designs, Psychon Bull Rev, 1(4), 476-490.
[2] Cousineau (2005), Confidence intervals in within-subject designs: a
    simpler solution to Loftus & Masson's method, Tutorials in Quantitative
    Psychology, 1(10, 42-45.
[3] Morey (2008). Intervals from Normalized Data: A correction to
    Cousineau (2005), Tutorials in Quantitative Psychology, 4(2), 61-64
[4] Franz & Loftus (2012). Standard errors and confidence intervals in
    within-subject designs: Generalizing Loftus & Masson (1994) and avoiding
    the biases of alternative accounts, Psychon Bull Rev, 19, 395-404
