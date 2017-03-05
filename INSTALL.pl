#!/usr/bin/env perl -w
use strict;

use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use File::Spec;
use Env qw(PATH);

BEGIN {

    $ENV{FP_HOME} = "$FindBin::Bin";

}

my $FPROOT = "$FindBin::RealBin";
my $all;
my $prep_genbank = "gunzip ".$FPROOT."/bin/GenBank_Plastomes.gz";
system($prep_genbank);
print "Thank you for dowloading the Fast-Plast pipeline. If you have not looked at the dependencies for Fast-Plast, please visit https://github.com/mrmckain/Fast-Plast and download them.\n\n";


print "Do you have all dependecies installed? Yes/No: ";

my $answer = <STDIN>;
chomp ($answer);
my $control_file =$FPROOT . "/fast-plast.pl";
if($answer =~ /n/i){

	print "\nDo you want me to try to install them? Yes/No/All: ";
	$answer = <STDIN>;
	chomp($answer);
	if($answer =~ /n/i){
		die "\nPlease install all dependencies prior to installing Fast-Plast.\n";
	}
	else{
		if($answer =~ /all/i){
				$all =1;
		}
		chdir("bin/");
		if($all){
			$answer="yes";
		}
		else{
		print "\nOK. I am only going to try to install the Linux binaries. Let's do this one-by-one. Would you like me to install Trimmomatic? Yes/No: ";
			$answer = <STDIN>;
			chomp($answer);
		}

		my $trimmomatic;
		if($answer =~ /y/i){
		
			system("wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip");
			
			if(-e "Trimmomatic-0.36.zip"){
				system("unzip Trimmomatic-0.36.zip");
				$trimmomatic = $FPROOT . "/bin/Trimmomatic-0.36/";
				$trimmomatic = glob ("$trimmomatic/*.jar");
			}
			else{
				print "\nUnable to install Trimmomatic.\n";
			}

		}
		else{
			print "Please provide the absolute path to the Trimmomatic directory. If Trimmomatic is in your PATH already, just type PATH: ";

			my $trimmomatic = <STDIN>;
			chomp ($trimmomatic);

			if($trimmomatic =~ /path/i){
				$trimmomatic= <"trimmomatic*.jar">;
				if(!$trimmomatic){

					my @path = split(/:/, $PATH);
					for my $pot (@path){
						if($pot =~ /trimmomatic/i){
							$trimmomatic = glob ("$pot/*.jar");
						}
					}
					if(!$trimmomatic){
						die "\nSorry. I cannot locate the Trimmomatic jar file in your path\. Please check again.\n";
					}
				}	

			}
			else{
				my $temp_trim = $trimmomatic;
				$trimmomatic = glob ($trimmomatic."*.jar");
				if(!$trimmomatic){
					die "\nSorry. I cannot locate the Trimmomatic jar file in $temp_trim\. Please check again.\n";
				}
			}

			print "\nTrimmomatic jar file located: $trimmomatic\n";
		}
		my $tempf = read_file($control_file);
		$tempf =~ s/my \$TRIMMOMATIC\;/my \$TRIMMOMATIC=\"$trimmomatic\"\;/;
		write_file($control_file, $tempf);
		

		if($all){
			$answer="yes";
		}
		else{
		print "\nWould you like me to install bowtie2? Yes/No: ";
		$answer = <STDIN>;
		chomp($answer);
	}
		my $bowtie2;
		if($answer =~ /y/i){
			system("wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip");
			if(-e "bowtie2-2.2.9-linux-x86_64.zip"){
				system("unzip bowtie2-2.2.9-linux-x86_64.zip");
				$bowtie2 = $FPROOT . "/bin/bowtie2-2.2.9";
				$bowtie2 = glob ("$bowtie2/bowtie2");
			}
			else{
				print "\nUnable to install bowtie2.\n";
			}
		}
		else{
			print "\nPlease provide the absolute path to the bowtie2 executable. If bowtie2 is in your PATH already, just type PATH: ";
			my $bowtie2 = <STDIN>;
			chomp ($bowtie2);
			if($bowtie2 =~ /path/i){
				$bowtie2=<"bowtie2">;
				if(!$bowtie2){
					my @path = split(/:/, $PATH);
					for my $pot (@path){
						if($pot =~ /$bowtie2/i){
							$bowtie2 = glob ("$pot/bowtie2");
						}
					}
					if(!$bowtie2){
						die "\nSorry. I cannot locate the bowtie2 executable file in your path\. Please check again.\n";
					}
				}

			}
			else{
				my $temp_bow = $bowtie2;
				$trimmomatic = glob ("$bowtie2/bowtie2");
				if(!$bowtie2){
					die "\nSorry. I cannot locate the bowtie2 executable file in $temp_bow\. Please check again.\n";
				}
			}

			print "\nbowtie2 executable located: $bowtie2\n";
		}
		$tempf = read_file($control_file);
		$tempf =~ s/my \$BOWTIE2\;/my \$BOWTIE2=\"$bowtie2\"\;/;
		write_file($control_file, $tempf);
		if($all){
			$answer="yes";
		}
		else{
			print "\nWould you like me to install SPAdes? Yes/No: ";
		$answer = <STDIN>;
		chomp($answer);
		}
		my $spades;
		if($answer =~ /y/i){	
			system("wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0-Linux.tar.gz");
			if(-e "SPAdes-3.9.0-Linux.tar.gz"){
				system("tar -xvzf SPAdes-3.9.0-Linux.tar.gz");
				$spades = $FPROOT . "/bin/SPAdes-3.9.0-Linux";
				$spades = glob ("$spades/bin/spades.py");
			}
			else{
				print "\nUnable to install SPAdes.\n";
			}
		}
		else{
			print "\nPlease provide the absolute path to the SPAdes program directory. If SPAdes is in your PATH already, just type PATH: ";

			my $spades = <STDIN>;
			chomp ($spades);

			if($spades =~ /path/i){
				$spades=<"spades.py">;
				if(!$spades){
					my @path = split(/:/, $PATH);
					for my $pot (@path){
						if($pot =~ /spades/i){
							$spades = glob ("$pot/bin/spades.py");
						}
					}
					if(!$spades){
						die "\nSorry. I cannot locate the SPAdes executable file in your path\. Please check again.\n";
					}
				}
			}
			else{
				my $temp_spades = $spades;
				$spades = glob ("$spades/spades.py");
				if(!$spades){
					$spades = glob ("$temp_spades/bin/spades.py");
				}
				if(!$spades){
					die "\nSorry. I cannot locate the SPAdes executable file in $temp_spades\. Please check again.\n";
				}
			}

			print "\nSPAdes executable located: $spades\n";
		}
		 $tempf = read_file($control_file);
		$tempf =~ s/my \$SPADES\;/my \$SPADES=\"$spades\"\;/;
		write_file($control_file, $tempf);

		if($all){
			$answer="yes";
		}
		else{
		print "\nWould you like me to install bowtie1? Yes/No: ";
		$answer = <STDIN>;
		chomp($answer);
	}
		my $bowtie1;
		if($answer =~ /y/i){
			system("wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip");
			if(-e "bowtie-1.1.2-linux-x86_64.zip"){
				system("unzip bowtie-1.1.2-linux-x86_64.zip");
				$bowtie1 = $FPROOT . "/bin/bowtie-1.1.2";
				$bowtie1 = glob ("$bowtie1/bowtie");
			}
			else{
				print "\nUnable to install bowtie1.\n";
			}
		}
		else{
			print "\nPlease provide the absolute path to the bowtie1 executable. If bowtie1 is in your PATH already, just type PATH: ";
			my $bowtie1 = <STDIN>;
			chomp ($bowtie1);
			if($bowtie1 =~ /path/i){
				$bowtie1=<"bowtie">;
				if(!$bowtie1){
					my @path = split(/:/, $PATH);
					for my $pot (@path){
						if($pot =~ /bowtie/i){
							$bowtie1 = glob ("$pot/bowtie");
						}
					}
					if(!$bowtie1){
						die "\nSorry. I cannot locate the bowtie1 executable file in your path\. Please check again.\n";
					}
				}
			}
			else{
				my $temp_bow1 = $bowtie1;
				$bowtie1 = glob ("$bowtie1/bowtie");
				if(!$bowtie1){
					die "\nSorry. I cannot locate the bowtie1 executable file in $temp_bow1\. Please check again.\n";
				}
			}

			print "\nbowtie1 executable located: $bowtie1\n";
		}
		$tempf = read_file($control_file);
		$tempf =~ s/my \$BOWTIE1\;/my \$BOWTIE1=\"$bowtie1\"\;/;
		write_file($control_file, $tempf);
		if($all){
			$answer="yes";
		}
		else{
		print "\nWould you like me to install SSPACE? Yes/No: ";
		$answer = <STDIN>;
		chomp($answer);
	}
		my $sspace;
		if($answer =~ /y/i){	
			system("wget https://github.com/nsoranzo/sspace_basic/archive/v2.1.1.zip");
			if(-e "v2.1.1.zip"){
				system("unzip v2.1.1.zip");
				$sspace = $FPROOT . "/bin/sspace_basic-2.1.1/";
				$sspace = glob ("$sspace/SSPACE_Basic.pl");
			}
			else{
				print "\nUnable to install SSPACE.\n";
			}
		}
		else{
			print "\nPlease provide the absolute path to the SSPACE program directory. If SSPACE is in your PATH already, just type PATH: ";

			my $sspace = <STDIN>;
			chomp ($sspace);

			if($sspace =~ /path/i){
				$sspace=<"SSPACE_Basic.pl">;
				if(!$sspace){
					my @path = split(/:/, $PATH);
					for my $pot (@path){
						if($pot =~ /sspace/i){
							$sspace = glob ("$pot/SSPACE_Basic.pl");
						}
					}
					if(!$sspace){
						die "\nSorry. I cannot locate the SSPACE executable file in your path\. Please check again.\n";
					}
				}
			}
			else{
				my $temp_sspace = $sspace;
				$sspace = glob ("$sspace/SSPACE_Basic.pl");
				if(!$sspace){
					$sspace = glob ("$temp_sspace/SSPACE_Basic.pl");
				}
				if(!$sspace){
					die "\nSorry. I cannot locate the SSPACE executable file in $temp_sspace\. Please check again.\n";
				}
			}

			print "\nSSPACE executable located: $sspace\n";
		}
		 $tempf = read_file($control_file);
		$tempf =~ s/my \$SSPACE\;/my \$SSPACE=\"$sspace\"\;/;
		write_file($control_file, $tempf);

		if($all){
			$answer="yes";
		}
		else{
		print "\nWould you like me to install BLAST+? Yes/No: ";
		$answer = <STDIN>;
		chomp($answer);
	}
		my $blastn;
		if($answer =~ /y/i){
			
			system("wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz");
			
			if(-e "ncbi-blast-2.6.0+-x64-linux.tar.gz"){
				system("tar -xvzf ncbi-blast-2.6.0+-x64-linux.tar.gz");
				$blastn = $FPROOT . "/bin/ncbi-blast-2.6.0+";
				$blastn = glob ("$blastn/bin/");
			}
			else{
				print "\nUnable to install BLAST+.\n";
			}
		}
		else{
			print "\nPlease provide the absolute path to the blastn executable. If blastn is in your PATH already, just type PATH: ";

			my $blastn = <STDIN>;
			chomp ($blastn);

			if($blastn =~ /path/i){
				$blastn=<"blastn">;
				if(!$blastn){
					my @path = split(/:/, $PATH);
					for my $pot (@path){
						if($pot =~ /blast/i){
							$blastn = glob ("$pot/blastn");
							if(!$blastn){
								$blastn = glob ("$pot/bin/");
							}
						}
					}
					if(!$blastn){
						die "\nSorry. I cannot locate the blastn executable in your path\. Please check again.\n";
					}
				}
			}
			else{
				my $temp_blastn = $blastn;
				$blastn = glob ("$blastn/blastn");
				if(!$blastn){
					$blastn = glob ("$temp_blastn/bin");
				}
				if(!$blastn){
					die "\nSorry. I cannot locate the blastn executable in $temp_blastn\. Please check again.\n";
				}
			}
		}
		
		 $tempf = read_file($control_file);
		$tempf =~ s/my \$BLAST\;/my \$BLAST=\"$blastn\"\;/;
		write_file($control_file, $tempf);
		
		if($all){
			$answer="yes";
		}
		else{
		print "\nWould you like me to install Jellyfish 2? Yes/No: ";
		$answer = <STDIN>;
		chomp($answer);
	}
		my $jellyfish;
		if($answer =~ /y/i){
			
			system("wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz");
			
			if(-e "jellyfish-2.2.6.tar.gz"){
				system("tar -xvzf jellyfish-2.2.6.tar.gz");
				$jellyfish = $FPROOT . "/bin/jellyfish-2.2.6/";
				chdir($jellyfish);
				system("./configure");
				system("make");
				$jellyfish = glob ($jellyfish."bin/jellyfish");
				chdir("../");
			}
			else{
				print "\nUnable to install Jellyfish2.\n";
			}
		}
		else{
			print "\nPlease provide the absolute path to the Jellyfish 2 executable. If Jellyfish is in your PATH already, just type PATH: ";

			my $jellyfish = <STDIN>;
			chomp ($jellyfish);

			if($jellyfish =~ /path/i){
				$jellyfish=<"jellyfish">;
				if(!$jellyfish){
					my @path = split(/:/, $PATH);
					for my $pot (@path){
						if($pot =~ /jellyfish/i){
							$jellyfish = glob ("$pot/jellyfish");
							if(!$jellyfish){
								$jellyfish = glob ("$pot/bin/jellyfish");
							}
						}
					}
					if(!$jellyfish){
						die "\nSorry. I cannot locate the jellyfish executable in your path\. Please check again.\n";
					}
				}
			}
			else{
				my $temp_jellyfish = $jellyfish;
				$jellyfish = glob ("$jellyfish/jellyfish");
				if(!$jellyfish){
					$jellyfish = glob ("$temp_jellyfish/bin/jellyfish");
					if(!$jellyfish){
						die "\nSorry. I cannot locate the jellyfish executable file in $temp_jellyfish\. Please check again.\n";
					}
				}
			}
			}
			print "\njellyfish executable located: $jellyfish\n";
			$tempf = read_file($control_file);
			$tempf =~ s/my \$JELLYFISH\;/my \$JELLYFISH=\"$jellyfish\"\;/;
			write_file($control_file, $tempf);

	}
	if($all){
			$answer="yes";
		}
	else{
	print "Would you like me to compile afin? Yes or No: ";

$answer = <STDIN>;
chomp($answer);
}
if($answer =~ /y/i){
	chdir("$FPROOT/afin");
	`make`;
	if(! -e "afin"){
		die "Could not compile afin.  You might be missing a required library.  Check that you have GCC 4.5+ and go to https://github.com/mrmckain/Fast-Plast/tree/master/afin for more information.\n";
	}
	print "\nAfin installed successfully.\n";
	chdir ("../");
}
print "Fast-Plast installation complete.  See manual for directions on how to get started.\n"
}
else{
	print "Great! We still need to set up Fast-Plast. Follow the promps.";




		print "To get started, we are going to set the paths to Trimmomatic, bowtie2, SPAdes, bowtie1, SSPACE, and BLAST+.  After set up for the main pipeline, we will try to set up the Coverage Analysis pipeline.\n\n";

		print "Please provide the absolute path to the Trimmomatic directory. If Trimmomatic is in your PATH already, just type PATH: ";

		my $trimmomatic = <STDIN>;
		chomp ($trimmomatic);

		if($trimmomatic =~ /path/i){
			$trimmomatic= <"trimmomatic*.jar">;
			
		if(!$trimmomatic){

			my @path = split(/:/, $PATH);
			for my $pot (@path){
			if($pot =~ /trimmomatic/i){
				$trimmomatic = glob ("$pot/*.jar");
				
			}
		}
		if(!$trimmomatic){
			die "\nSorry. I cannot locate the Trimmomatic jar file in your path\. Please check again.\n";
		}
		
	}	
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$TRIMMOMATIC\;/my \$TRMIMMOMATIC=\"$trimmomatic\"\;/;
		write_file($control_file, $tempf);

}
else{
	my $temp_trim = $trimmomatic;
	$trimmomatic = glob ("$trimmomatic/*.jar");
	if(!$trimmomatic){
		die "\nSorry. I cannot locate the Trimmomatic jar file in $temp_trim\. Please check again.\n";
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$TRIMMOMATIC\;/my \$TRMIMMOMATIC=\"$trimmomatic\"\;/;
		write_file($control_file, $tempf);
}

print "\nTrimmomatic jar file located: $trimmomatic\n";

print "\nPlease provide the absolute path to the bowtie2 executable. If bowtie2 is in your PATH already, just type PATH: ";

my $bowtie2 = <STDIN>;
chomp ($bowtie2);

if($bowtie2 =~ /path/i){
	$bowtie2=<"bowtie2">;
	if(!$bowtie2){
		my @path = split(/:/, $PATH);
		for my $pot (@path){
			if($pot =~ /$bowtie2/i){
				$bowtie2 = glob ("$pot/bowtie2");
			}
		}
		if(!$bowtie2){
			die "\nSorry. I cannot locate the bowtie2 executable file in your path\. Please check again.\n";
		}
		
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$BOWTIE2\;/my \$BOWTIE2=\"$bowtie2\"\;/;
		write_file($control_file, $tempf);

}
else{
	my $temp_bow = $bowtie2;
	$trimmomatic = glob ("$bowtie2/bowtie2");
	if(!$bowtie2){
		die "\nSorry. I cannot locate the bowtie2 executable file in $temp_bow\. Please check again.\n";
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$BOWTIE2\;/my \$BOWTIE2=\"$bowtie2\"\;/;
		write_file($control_file, $tempf);
}

print "\nbowtie2 executable located: $bowtie2\n";


print "\nPlease provide the absolute path to the SPAdes program directory. If SPAdes is in your PATH already, just type PATH: ";

my $spades = <STDIN>;
chomp ($spades);

if($spades =~ /path/i){
	$spades=<"spades.py">;
	if(!$spades){
		my @path = split(/:/, $PATH);
		for my $pot (@path){
			if($pot =~ /spades/i){
				$spades = glob ("$pot/bin/spades.py");
			}
		}
		if(!$spades){
			die "\nSorry. I cannot locate the SPAdes executable file in your path\. Please check again.\n";
		}
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$SPADES\;/my \$SPADES=\"$spades\"\;/;
		write_file($control_file, $tempf);
}
else{
	my $temp_spades = $spades;
	$spades = glob ("$spades/spades.py");
	if(!$spades){
		$spades = glob ("$temp_spades/bin/spades.py");
	}
	if(!$spades){
		die "\nSorry. I cannot locate the SPAdes executable file in $temp_spades\. Please check again.\n";
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$SPADES\;/my \$SPADES=\"$spades\"\;/;
		write_file($control_file, $tempf);
}

print "\nSPAdes executable located: $spades\n";

print "\nPlease provide the absolute path to the bowtie1 executable. If bowtie1 is in your PATH already, just type PATH: ";

my $bowtie1 = <STDIN>;
chomp ($bowtie1);

if($bowtie1 =~ /path/i){
	$bowtie1=<"bowtie">;
	if(!$bowtie1){
		my @path = split(/:/, $PATH);
		for my $pot (@path){
			if($pot =~ /bowtie/i){
				$bowtie1 = glob ("$pot/bowtie");
			}
		}
		if(!$bowtie1){
			die "\nSorry. I cannot locate the bowtie1 executable file in your path\. Please check again.\n";
		}
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$BOWTIE1\;/my \$BOWTIE1=\"$bowtie1\"\;/;
		write_file($control_file, $tempf);
}
else{
	my $temp_bowtie1 = $bowtie1;
	$bowtie1 = glob ("$bowtie1/bowtie");
	if(!$bowtie1){
		$bowtie1 = glob ("$temp_bowtie1/bowtie");
	}
	if(!$bowtie1){
		die "\nSorry. I cannot locate the bowtie1 executable file in $temp_bowtie1\. Please check again.\n";
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$BOWTIE1\;/my \$BOWTIE1=\"$bowtie1\"\;/;
		write_file($control_file, $tempf);
}
print "\nBowtie1 executable located: $bowtie1\n";


print "\nPlease provide the absolute path to the SSPACE executable. If SSPACE is in your PATH already, just type PATH: ";

my $sspace = <STDIN>;
chomp ($sspace);

if($sspace =~ /path/i){
	$sspace=<"SSPACE_Basic.pl">;
	if(!$sspace){
		my @path = split(/:/, $PATH);
		for my $pot (@path){
			if($pot =~ /SSPACE/i){
				$sspace = glob ("$pot/SSPACE_Basic.pl");
			}
		}
		if(!$sspace){
			die "\nSorry. I cannot locate the SSPACE executable file in your path\. Please check again.\n";
		}
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$SSPACE\;/my \$SSPACE=\"$sspace\"\;/;
		write_file($control_file, $tempf);
}
else{
	my $temp_sspace = $sspace;
	$sspace = glob ("$sspace/SSPACE_Basic.pl");
	if(!$sspace){
		$sspace = glob ("$temp_sspace/SSPACE_Basic.pl");
	}
	if(!$sspace){
		die "\nSorry. I cannot locate the SSPACE executable file in $temp_sspace\. Please check again.\n";
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$SSPACE\;/my \$SSPACE=\"$sspace\"\;/;
		write_file($control_file, $tempf);
}

print "\nSSPACE executable located: $sspace\n";

print "\nPlease provide the absolute path to the blastn executable. If blastn is in your PATH already, just type PATH: ";

my $blastn = <STDIN>;
chomp ($blastn);

if($blastn =~ /path/i){
	$blastn=<"blastn">;
	if(!$blastn){
		my @path = split(/:/, $PATH);
		for my $pot (@path){
			if($pot =~ /blast/i){
				$blastn = glob ("$pot/blastn");
				if(!$blastn){
					$blastn = glob ("$pot/bin/blastn");
				}
			}
		}
		if(!$blastn){
			die "\nSorry. I cannot locate the blastn executable in your path\. Please check again.\n";
		}
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$BLAST\;/my \$BLAST=\"$blastn\"\;/;
		write_file($control_file, $tempf);
}
else{
	my $temp_blastn = $blastn;
	$blastn = glob ("$blastn/blastn");
	if(!$blastn){
		$blastn = glob ("$temp_blastn/bin/blastn");
	}
	if(!$blastn){
		die "\nSorry. I cannot locate the blastn executable in $temp_blastn\. Please check again.\n";
	}
	my $tempf = read_file($control_file);
		$tempf =~ s/my \$BLAST\;/my \$BLAST=\"$blastn\"\;/;
		write_file($control_file, $tempf);
}

print "\nblastn executable located: $blastn\n";


print "Would you like me to compile afin? Yes or No: ";

$answer = <STDIN>;
chomp($answer);

if($answer =~ /y/i){
	chdir("../afin");
	`make`;
	if(! -e "afin"){
		die "Could not compile afin.  You might be missing a required library.  Check that you have GCC 4.5+ and go to https://github.com/mrmckain/Fast-Plast/tree/master/afin for more information.\n";
	}
	print "\nAfin installed successfully.\n";
	chdir ("../bin");
}



`perl -pi -e "s/my \$TRIMMOMATIC;/my \$TRIMMOMATIC = $trimmomatic\;/" $FPROOT/fast-plast.pl`;
`perl -pi -e "s/my \$BOWTIE2;/my \$BOWTIE2 = $bowtie2\;/" $FPROOT/fast-plast.pl`;
`perl -pi -e "s/my \$SPADES;/my \$SPADES = $spades\;/" $FPROOT/fast-plast.pl`;
`perl -pi -e "s/my \$BLAST;/my \$BLAST = $blastn\;/" $FPROOT/fast-plast.pl`;




	print "\nPlease provide the absolute path to the Jellyfish 2 executable. If Jellyfish is in your PATH already, just type PATH: ";

	my $jellyfish = <STDIN>;
	chomp ($jellyfish);

	if($jellyfish =~ /path/i){
		$jellyfish=<"jellyfish">;
		if(!$jellyfish){
			my @path = split(/:/, $PATH);
			for my $pot (@path){
				if($pot =~ /jellyfish/i){
					$jellyfish = glob ("$pot/jellyfish");
					if(!$jellyfish){
						$jellyfish = glob ("$pot/bin/jellyfish");
					}
				}
			}
			if(!$jellyfish){
				die "\nSorry. I cannot locate the jellyfish executable in your path\. Please check again.\n";
			}
		}
		print "\njellyfish executable located: $jellyfish\n";
			my $tempf = read_file($control_file);
			$tempf =~ s/my \$JELLYFISH\;/my \$JELLYFISH=\"$jellyfish\"\;/;
			write_file($control_file, $tempf);
	}
	else{
		my $temp_jellyfish = $jellyfish;
		$jellyfish = glob ("$jellyfish/jellyfish");
		if(!$jellyfish){
			$jellyfish = glob ("$temp_jellyfish/bin/jellyfish");
			if(!$jellyfish){
				die "\nSorry. I cannot locate the jellyfish executable file in $temp_jellyfish\. Please check again.\n";
			}
		}
		print "\njellyfish executable located: $jellyfish\n";
			my $tempf = read_file($control_file);
			$tempf =~ s/my \$JELLYFISH\;/my \$JELLYFISH=\"$jellyfish\"\;/;
			write_file($control_file, $tempf);
	}

	print "\njellyfish executable located: $jellyfish\n";

	`perl -pi -e "s/my \$JELLYFISH;/my \$JELLYFISH = $jellyfish\;/" $FPROOT/Coverage_Analysis/coverage.pl`;

print "Fast-Plast installation complete.  See manual for directions on how to get started.\n"
}


sub read_file {
    my ($filename) = @_;
 
    open my $in, '<:encoding(UTF-8)', $filename or die "Could not open '$filename' for reading $!";
    local $/ = undef;
    my $all = <$in>;
    close $in;
 
    return $all;
}

sub write_file {
    my ($filename, $content) = @_;
 
    open my $out, '>:encoding(UTF-8)', $filename or die "Could not open '$filename' for writing $!";;
    print $out $content;
    close $out;
 
    return;
}


