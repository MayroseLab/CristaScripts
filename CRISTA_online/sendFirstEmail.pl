#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;

my $toEmail;
my $jobId;

GetOptions(
	"toEmail=s"        => \$toEmail,
	"id=s"     => \$jobId
);

my $from 		= 'evolseq@post.tau.ac.il';
my $subject	= "Your CRISTA job #$jobId is being processed";

my $resultsLink = 'http://crista.tau.ac.il/results/';
$resultsLink = $resultsLink."$jobId";
$resultsLink = $resultsLink."/output.php";

my $message   = "Dear CRISTA user,\n\n";
$message  .= "Your CRISTA job #$jobId is being processed.\n";
$message  .= "We will email you again when processing is done.\n\n";
$message  .= "View your job\'s progress and results here:\n";
$message  .= "$resultsLink\n\n";

$message  .= "Thanks you for using CRISTA!\n";
$message  .= "CRISTA team.\n";

open(MAIL, "|/usr/sbin/sendmail -t");
 
# Email Header
print MAIL "To: $toEmail\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject\n\n";

# Email Body
print MAIL $message;

close(MAIL);
print "Email Sent Successfully to: $toEmail\n";

