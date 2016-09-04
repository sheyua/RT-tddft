# LaTeX2HTML 2008 (1.71)
# Associate labels original text with physical files.


1;


# LaTeX2HTML 2008 (1.71)
# labels from external_latex_labels array.


$key = q//;
$external_latex_labels{$key} = q|\fi|; 
$noresave{$key} = "$nosave";

$key = q/_newlabelxx/;
$external_latex_labels{$key} = q|\ifx|; 
$noresave{$key} = "$nosave";

1;

