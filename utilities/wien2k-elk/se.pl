#!/usr/bin/perl
if (@ARGV<1)
{
print "Script for conversion Wien2K struct files to spacegroup.in files used in exciting/elk\n";
print "The only argument is the name of the struct file, spacegroup.in is written to standard output\n";
print "Use at your own risk.  \n Jerzy Goraus 2009\n";
exit;
};
open(STR,$ARGV[0]) or die "can't open ",$ARGV[0]," file\n ";
@STRUCT=<STR>;
($TMP,$GROUP)=split(/\_/,$STRUCT[1]);
($GROUP,$TMP)=split(/\n/,$GROUP);
($GROUP,$TMP)=split(/\ /,$GROUP);
($TMP,$LA,$LB,$LC,$AA,$AB,$AG)=split(/\s+/,$STRUCT[3]);
$ATCNT=0;
print " '",$GROUP,"\'\n ",$LA," ",$LB," ",$LC,"\n ",$AA," ",$AB," ",$AG,"\n 1  1  1\n .false.\n ";
 
foreach (@STRUCT)
{
$BUF=$_;
if (/ATOM/) 
 { 
  
  ($TMP,$X,$Y,$Z)=split(/=/,$BUF);
  ($X,$TMP)=split(/\ /,$X);
  ($Y,$TMP)=split(/\ /,$Y);
  $ATCNT+=1;
#  print "\n",$X, " ", $Y," ", $Z;
  $RESULT = sprintf(" %lf %lf %lf\n", $X,$Y,$Z);
 };
if (/NPT/) 
 {
 $ATL=substr($_,0,2);
 $_=$ATL;
 s/\s//g;
 $ATL=$_;
# print $ATL;
 $TAB{$ATL}=$TAB{$ATL}.$RESULT;
 $TABN{$ATL}+=1;
 };
};
$NUM=keys(TAB);
print " ",$NUM ,"\n";
foreach (keys(TAB))
{
 print " '",$_,"' '",$_,".in'\n";
 print " ",$TABN{$_},"\n";
 print $TAB{$_};
}


