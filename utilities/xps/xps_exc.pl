#!/usr/bin/perl 
print "Script for calculating VB-XPS spectra from  PDOS* files (with option dosmsum and lmi rep switched off) and crossection \n";
print "for photoemission read from file CROSEC (containing s,p,d,f crossections). CROSEC should be placed in the same directory\n ";
print "(you can change this in this script at line with comment - Here enter CROSEC). Some elements I used are already in CROSEC\n ";
print "(for E=1486eV, Al Kalpha) but in principle depending on the energy of XPS measurement you should use  your own CROSEC,\n";
print " take the needed values e.q. from  J. J. Yeh and I. Lindau, Atomic Data and Nuclear Data Tables 32 1 (1985).\n";
print "The PDOSES multiplied by crossection and summed are written to _CTDOS file. In order to compare the result with experimental VB-XPS spectra\n" ;
print  "one should convolute the _cTDOS  with Lorentzian and Fermi function. For example with the use of conv2.c (also provided)\n";
print "The script also produces s,p,d,f sums of PDOSes summed across m number (files ending with four digit number) and \n";
print "summed across all Wyckoff positions for a considered element (files _elname.dat). In this script line starting with FACTOR\n" ;
print "can be used to convert between Ha units and eV, for convolution purposes it should be left as is, so the script\n";
print "reads the PDOSes and outputs the abscissa in units of eV.\n";
print "The script produces also .gnu files which  can be plotted with gnuplot by running gnuplot *.gnu .\n";
print "I provide this script for convenience of use. I take no responsibility for wrong results, damages etc. \n";
print "so you use it at your own risk. Of course you can modify this script according to your needs.\n";
print " Jerzy Goraus, (2009)\n"; 
print "----------------------------------------------------------------------------------------------------------------------\n\n";

# corrected bug in line 127..129, 20 Oct 2010
$GR=0;
$TERM="set term png med";
$SUF="png";
$PREC="-0.50 0.002 1.500 0.003";
$OPT="set key left";
$SPOL=0;
$FACTOR=27.211;
if ($ARGV[1] eq "-g") { $GR=1;}
$atidx=1;
open(STR,  "elk.in")   or die "Can't open elk.in\n";
@str=<STR>;
close(STR);
$atcnt=0;
print "\nElements : ";
for  ($i=0; $i<$#str; $i++)
{
$_=$str[$i];
if (/atoms\n/)
{
$i++;
$_=$str[$i];
s/\s//g;
($ATCNT,$TMP)=split(/:/);
}

if (/spinpol/)
{
$i++;
$_=$str[$i];
if (/true/) { $SPOL=1;}
}

if (/in\'/)
{
($ELNAME,$TMP)=split(/\./);
$_=$ELNAME;
s/\'//;
$ELNAME=$_;
#print "\n",$ELNAME;
$_=$str[$i+1];
s/\s//g;
($MULT,$TMP)=split(/:/);
push @atoms, $ELNAME;
print $ELNAME," (",$MULT,") ";

$mtab{$ELNAME}=$MULT;
$idxtab{$ELNAME}=$atidx;
$atidx++;
#print "|",$MULT,"|";

}
}

if ($SPOL ) { print "\nMagnetic case\n"; } else { print "\nNon-magnetic case\n";}
#print $ATCNT;
open(CRS,"./CROSEC")   or die "Can't open CROSEC";  # Here enter the position of CROSEC
print "\nRead CROSEC";
@crs=<CRS>;
close(CRS);
foreach ($i=1 ; $i<@crs; $i++)
{
$_=$crs[$i];
if (/ATOM=/)
{
$A=substr($_,5);
$A=~s/\s//g;
} else
{
$CR{$A}=$_;
}
}

print "\n";
 $NDOS=0;
foreach (@atoms)
{

 $ELNAME=$_;
 for ($i=1; $i<=$mtab{$_}; $i++)
 {
   $fname=sprintf("PDOS_S%02d_A%04d.OUT", $idxtab{$_}, $i);
   print "\nRead : ",$fname  ;
   open(pdos,  $fname)   or die "Can't open $fname\n";
   @pdos=<pdos>;
   close(pdos);
   if ($NDOS==0) 
   {
   foreach (@pdos) 
   {
    if (/\d/) {$NDOS++}; 
   }
   $NDOS/=32;
   $NDOS++;
   #- init to 0
      for ($j=0; $j<$NDOS-1; $j++) 
	{
    	 $atw_su[$j]=0; $atw_pu[$j]=0; $atw_du[$j]=0; $atw_fu[$j]=0; 
    	 $atw_sd[$j]=0; $atw_pd[$j]=0; $atw_dd[$j]=0; $atw_fd[$j]=0; 
    		 
    	 $att_su[$j]=0; $att_pu[$j]=0; $att_du[$j]=0; $att_fu[$j]=0; 
    	 $att_sd[$j]=0; $att_pd[$j]=0; $att_dd[$j]=0; $att_fd[$j]=0; 
    	 
         $attc_su[$j]=0; $attc_pu[$j]=0; $attc_du[$j]=0; $attc_fu[$j]=0; 
    	 $attc_sd[$j]=0; $attc_pd[$j]=0; $attc_dd[$j]=0; $attc_fd[$j]=0; 
    	 
        
        }

   } 
   $atfout=sprintf("_%s_%04d.dat",$ELNAME,$i);
   open(atpart,  ">".$atfout)   or die "Can't open output file $atfout\n";
   if ($SPOL ) { print atpart "# spin up - s,p,d,f and spin down s,p,d,f\n"; } else { print atpart "# s, p, d, f\n";}
   
   
   for ($q=0; $q<$#pdos; $q++) 
   { 
   ($TMP,$pdosx[$q],$pdosy[$q])=split(/\s+/,$pdos[$q]);
   $pdosx[$q]*=$FACTOR; $pdosy[$q]/=$FACTOR;
   }

   for ($j=0; $j<$NDOS-1; $j++)
   {
    $at_su[$j]=$pdosy[$j];
    if ($SPOL ) {$at_sd[$j]=$pdosy[$j+16*$NDOS];}
    $at_pu[$j]=0; $at_pd[$j]=0; 
    $at_du[$j]=0; $at_dd[$j]=0; 
    $at_fu[$j]=0; $at_fd[$j]=0; 
    for ($q=1; $q<4; $q++) 
    {
        $at_pu[$j]+=$pdosy[$j+$q*$NDOS];
        if ($SPOL ) {$at_pd[$j]+=$pdosy[$j+16*$NDOS+$q*$NDOS];}
    }
    for ($q=4; $q<9; $q++) 
    {
        $at_du[$j]+=$pdosy[$j+$q*$NDOS];
        if ($SPOL ) {$at_dd[$j]+=$pdosy[$j+16*$NDOS+$q*$NDOS];}
    }
    for ($q=9; $q<16; $q++) 
    {
        $at_fu[$j]+=$pdosy[$j+$q*$NDOS];
        if ($SPOL ) {$at_fd[$j]+=$pdosy[$j+16*$NDOS+$q*$NDOS];}
    }
    

    print atpart $pdosx[$j], " ",$at_su[$j]," ",$at_pu[$j]," ",$at_du[$j]," ",$at_fu[$j];
    if ($SPOL ) { print atpart " ",$at_sd[$j]," ",$at_pd[$j]," ",$at_dd[$j]," ",$at_fd[$j],"\n";} else {print atpart "\n";}
    
    $atw_su[$j]+=$at_su[$j]; $atw_pu[$j]+=$at_pu[$j]; $atw_du[$j]+=$at_du[$j]; $atw_fu[$j]+=$at_fu[$j];
    $atw_sd[$j]+=$at_sd[$j]; $atw_pd[$j]+=$at_pd[$j]; $atw_dd[$j]+=$at_dd[$j]; $atw_fd[$j]+=$at_fd[$j];
   }
   close atpart;
   print "\nWritten ",$atfout;
   #--
        $gatfout=sprintf("_g_%s_%04d",$ELNAME,$i);
   	open(GNU,">".$gatfout.".gnu") or die "can't open $gatfout\n";
	print GNU $TERM,"\n",$OPT,"\n";
	print GNU "set output \"$gatfout.$SUF\"\n";
	print GNU "set xlabel \"Eb [eV]\"\n";
	print GNU "set ylabel \"DOS [st/eV f.u.]\"\n";
	print GNU "set title  \"$ELNAME - $i\"\n";
	print GNU "plot \"$atfout\" using 1:2 smooth uniq t 's' ";
	printf GNU ", \"$atfout\" using 1:3 smooth uniq t \'p\' ";
	printf GNU ", \"$atfout\" using 1:4 smooth uniq t \'d\' ";
	printf GNU ", \"$atfout\" using 1:5 smooth uniq t \'f\' ";
	if ($SPOL ) 
	        { 
	        print GNU ", \"$atfout\" using 1:6 smooth uniq t \'s\' ";
	        print GNU ", \"$atfout\" using 1:7 smooth uniq t \'p\' ";
	        print GNU ", \"$atfout\" using 1:8 smooth uniq t \'d\' ";
	        print GNU ", \"$atfout\" using 1:9 smooth uniq t \'f\' ";
	        }
	close(GNU);
        print "\nWritten ",$gatfout,".gnu";	
   #--

   
   
 }
   $atwfout=sprintf("_%s.dat",$ELNAME);
#at - atomic DOS, atw - given specie dos, att - total dos , attc - total dos with crossection
   open(atwpart,  ">".$atwfout)   or die "Can't open output file $atwfout \n";
      if ($SPOL ) { print atwpart "# spin up - s,p,d,f and spin down s,p,d,f\n"; } else { print atwpart "# s, p, d, f\n";}
   if ($CR{$ELNAME} ne "")
   {  ($Ks,$Kp,$Kd,$Kf)=split(/\,/,$CR{$ELNAME}); } else {$Ks=$Kp=$Kd=$Kf=0; print "\n----------------- No $ELNAME in CROSEC, Using 0 as crossection\n";};
 
   
   for ($j=0; $j<$NDOS-1; $j++)
   {
   $att_su[$j]+=$atw_su[$j]; $att_pu[$j]+=$atw_pu[$j]; $att_du[$j]+=$atw_du[$j]; $att_fu[$j]+=$atw_fu[$j];
   $att_sd[$j]+=$atw_sd[$j]; $att_pd[$j]+=$atw_pd[$j]; $att_dd[$j]+=$atw_dd[$j]; $att_fd[$j]+=$atw_fd[$j];
   

   $attc_su[$j]+=$atw_su[$j]*$Ks; $attc_pu[$j]+=$atw_pu[$j]*$Kp; $attc_du[$j]+=$atw_du[$j]*$Kd; $attc_fu[$j]+=$atw_fu[$j]*$Kf;
   $attc_sd[$j]+=$atw_sd[$j]*$Ks; $attc_pd[$j]+=$atw_pd[$j]*$Kp; $attc_dd[$j]+=$atw_dd[$j]*$Kd; $attc_fd[$j]+=$atw_fd[$j]*$Kf;
   print atwpart $pdosx[$j], " ",$atw_su[$j]," ",$atw_pu[$j]," ",$atw_du[$j]," ",$atw_fu[$j];
   if ($SPOL ) { print atwpart " ",$atw_sd[$j]," ",$atw_pd[$j]," ",$atw_dd[$j]," ",$atw_fd[$j],"\n";} else {print atwpart "\n";}
   }
   close atwpart;
   #-
        $gatwfout=sprintf("_g_%s",$ELNAME);
   	open(GNU,">".$gatwfout.".gnu") or die "can't open $gatwfout\n";
	print GNU $TERM,"\n",$OPT,"\n";
	print GNU "set output \"$gatwfout.$SUF\"\n";
	print GNU "set xlabel \"Eb [eV]\"\n";
	print GNU "set ylabel \"DOS [st/eV f.u.]\"\n";
	print GNU "set title  \"$ELNAME \"\n";
	print GNU "plot \"$atwfout\" using 1:2 smooth uniq t 's' ";
	printf GNU ", \"$atwfout\" using 1:3 smooth uniq t \'p\' ";
	printf GNU ", \"$atwfout\" using 1:4 smooth uniq t \'d\' ";
	printf GNU ", \"$atwfout\" using 1:5 smooth uniq t \'f\' ";
	if ($SPOL ) 
	        { 
	        print GNU ", \"$atwfout\" using 1:6 smooth uniq t \'s\' ";
	        print GNU ", \"$atwfout\" using 1:7 smooth uniq t \'p\' ";
	        print GNU ", \"$atwfout\" using 1:8 smooth uniq t \'d\' ";
	        print GNU ", \"$atwfout\" using 1:9 smooth uniq t \'f\' ";
	        }

	close(GNU);
        print "\nWritten ",$gatwfout,".gnu";	
   
   
   #-
   print "\nWritten ",$atwfout;
   
   for ($j=0; $j<$NDOS-1; $j++) 
   {
        $atw_su[$j]=0; $atw_pu[$j]=0; $atw_du[$j]=0; $atw_fu[$j]=0; 
        $atw_sd[$j]=0; $atw_pd[$j]=0; $atw_dd[$j]=0; $atw_fd[$j]=0; 
   }
 
 
}

   open(att,  ">_TDOS.dat")   or die "Can't open output file TDOS.dat \n";
   open(attc,  ">_cTDOS.dat")   or die "Can't open output file cTDOS.dat \n";
   for ($j=0; $j<$NDOS-1; $j++)
   {
   $u=$att_su[$j]+$att_pu[$j]+$att_du[$j]+$att_fu[$j];
   $d=$att_sd[$j]+$att_pd[$j]+$att_dd[$j]+$att_fd[$j];
   $uc=$attc_su[$j]+$attc_pu[$j]+$attc_du[$j]+$attc_fu[$j];
   $dc=$attc_sd[$j]+$attc_pd[$j]+$attc_dd[$j]+$attc_fd[$j];
   if ($SPOL ) 
   {
        print att $pdosx[$j], " ",$u-$d," ",$u," ",$d,"\n";
        print attc $pdosx[$j], " ",$uc-$dc," ",$uc," ",$dc,"\n";
   } else
   {
        print att $pdosx[$j], " ",$u,"\n";
        print attc $pdosx[$j], " ",$u,"\n";
   } 
   }
   close att;
   print "\nWritten _TDOS.dat";   
   close attc;
   print "\nWritten _cTDOS.dat";   
   print "\n";
   
